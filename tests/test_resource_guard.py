import argparse
from pathlib import Path
from types import SimpleNamespace

import pytest

from scripts import resource_guard


def test_run_argument_parsing_preserves_command():
    args = resource_guard.build_parser().parse_args([
        "run", "--path", "/data", "--min-free-gb", "35", "--memory-gb", "16",
        "--", "python", "analysis.py", "--sample", "A",
    ])

    assert args.path == Path("/data")
    assert args.min_free_gb == 35
    assert args.memory_gb == 16
    assert resource_guard.normalize_command(args.command) == [
        "python", "analysis.py", "--sample", "A",
    ]


def test_download_argument_parsing():
    args = resource_guard.build_parser().parse_args([
        "download", "--url", "https://example.org/file.h5ad", "--output", "/data/file.h5ad",
        "--min-free-gb", "35", "--max-file-gb", "10",
    ])

    assert args.output == Path("/data/file.h5ad")
    assert args.min_free_gb == 35
    assert args.max_file_gb == 10


@pytest.mark.parametrize("value", ["0", "-1", "nan", "inf", "not-a-number"])
def test_gb_parser_rejects_invalid_values(value):
    with pytest.raises(argparse.ArgumentTypeError):
        resource_guard.positive_gb(value)


@pytest.mark.parametrize("value,expected", [("0", 0), ("123", 123), (" 456 ", 456)])
def test_content_length_parsing(value, expected):
    assert resource_guard.parse_content_length(value) == expected


@pytest.mark.parametrize("value", [None, "", "unknown", "1.5", "-1"])
def test_content_length_parsing_fails_closed(value):
    with pytest.raises(resource_guard.ResourceGuardError, match="Content-Length"):
        resource_guard.parse_content_length(value)


def test_download_decision_uses_only_remaining_bytes():
    remaining = resource_guard.assess_download(
        content_length=8 * resource_guard.GIB,
        partial_size=4 * resource_guard.GIB,
        free_bytes=15 * resource_guard.GIB,
        min_free_gb=10,
        max_file_gb=10,
    )

    assert remaining == 4 * resource_guard.GIB


def test_download_decision_allows_exact_disk_reserve_and_size_limit():
    remaining = resource_guard.assess_download(
        content_length=5 * resource_guard.GIB,
        partial_size=resource_guard.GIB,
        free_bytes=14 * resource_guard.GIB,
        min_free_gb=10,
        max_file_gb=5,
    )

    assert remaining == 4 * resource_guard.GIB


def test_download_decision_rejects_excessive_file_size():
    with pytest.raises(resource_guard.ResourceGuardError, match="maximum file size"):
        resource_guard.assess_download(11 * resource_guard.GIB, 0, 100 * resource_guard.GIB, 10, 10)


def test_download_decision_rejects_disk_reserve_violation():
    with pytest.raises(resource_guard.ResourceGuardError, match="reserve required"):
        resource_guard.assess_download(
            8 * resource_guard.GIB, 0, 17 * resource_guard.GIB, 10, 10
        )


def test_download_decision_rejects_oversized_partial_file():
    with pytest.raises(resource_guard.ResourceGuardError, match="refusing to resume"):
        resource_guard.assess_download(
            4 * resource_guard.GIB, 5 * resource_guard.GIB, 100 * resource_guard.GIB, 10, 10
        )


def test_run_disk_check_fails_before_launch_when_below_minimum(monkeypatch):
    monkeypatch.setattr(
        resource_guard.shutil,
        "disk_usage",
        lambda path: SimpleNamespace(free=9 * resource_guard.GIB),
    )

    with pytest.raises(resource_guard.ResourceGuardError, match="insufficient free space"):
        resource_guard.check_free_space("/data", 10)


def test_systemd_command_construction_disables_swap():
    command = resource_guard.build_systemd_command(16.0, ["python", "analysis.py", "--fast"])

    assert command == [
        "systemd-run", "--user", "--scope", "--property=MemoryMax=16G",
        "--property=MemorySwapMax=0", "--", "python", "analysis.py", "--fast",
    ]


def test_curl_command_construction_enables_resume_fail_retry_and_size_cap(tmp_path):
    output = tmp_path / "partial.bin"
    command = resource_guard.build_curl_command(
        "https://example.org/file.bin", output, 10 * resource_guard.GIB
    )

    assert command[:2] == ["curl", "--fail"]
    assert command[command.index("--continue-at") + 1] == "-"
    assert command[command.index("--retry") + 1] == "5"
    assert "--retry-all-errors" in command
    assert command[command.index("--max-filesize") + 1] == str(10 * resource_guard.GIB)
    assert command[command.index("--output") + 1] == str(output)
    assert command[-1] == "https://example.org/file.bin"
