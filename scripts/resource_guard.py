#!/usr/bin/env python3
"""Guard heavy commands and resumable downloads with resource limits."""

import argparse
import contextlib
import fcntl
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import urllib.error
import urllib.request

GIB = 1024 ** 3
DEFAULT_LOCK_PATH = (
    Path(os.environ.get("XDG_RUNTIME_DIR", tempfile.gettempdir()))
    / "agescape-resource-guard-heavy-job.lock"
)


class ResourceGuardError(RuntimeError):
    pass


def positive_gb(value):
    try:
        size = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"invalid GB value: {value!r}") from exc
    if not math.isfinite(size) or size <= 0:
        raise argparse.ArgumentTypeError("GB value must be a finite number greater than zero")
    return size


def gb_to_bytes(size_gb):
    return int(size_gb * GIB)


def format_bytes(size):
    return f"{size / GIB:.2f} GB"


def format_gb(size_gb):
    return str(int(size_gb)) if size_gb.is_integer() else str(size_gb)


def normalize_command(command):
    if command and command[0] == "--":
        command = command[1:]
    if not command:
        raise ResourceGuardError("a command is required after --")
    return command


def build_systemd_command(memory_gb, command):
    command = normalize_command(command)
    return [
        "systemd-run",
        "--user",
        "--scope",
        f"--property=MemoryMax={format_gb(memory_gb)}G",
        "--property=MemorySwapMax=0",
        "--",
        *command,
    ]


def build_curl_command(
    url,
    output,
    max_file_bytes,
    connect_timeout=30,
    max_time=None,
    speed_limit=10240,
    speed_time=120,
):
    command = [
        "curl",
        "--fail",
        "--location",
        "--continue-at",
        "-",
        "--retry",
        "5",
        "--retry-delay",
        "5",
        "--retry-all-errors",
        "--connect-timeout",
        str(int(connect_timeout)),
        "--max-filesize",
        str(max_file_bytes),
    ]
    if max_time is not None and int(max_time) > 0:
        command.extend(["--max-time", str(int(max_time))])
    if speed_limit is not None and int(speed_limit) > 0 and speed_time is not None and int(speed_time) > 0:
        command.extend(["--speed-limit", str(int(speed_limit)), "--speed-time", str(int(speed_time))])
    command.extend([
        "--output",
        str(output),
        url,
    ])
    return command


def check_free_space(path, min_free_gb):
    try:
        free_bytes = shutil.disk_usage(path).free
    except OSError as exc:
        raise ResourceGuardError(f"cannot check free space at {path}: {exc}") from exc
    required_bytes = gb_to_bytes(min_free_gb)
    if free_bytes < required_bytes:
        raise ResourceGuardError(
            f"insufficient free space at {path}: {format_bytes(free_bytes)} available, "
            f"{format_bytes(required_bytes)} required"
        )
    return free_bytes


def parse_content_length(value):
    if value is None:
        raise ResourceGuardError("server response has no Content-Length; refusing download")
    try:
        content_length = int(value)
    except (TypeError, ValueError) as exc:
        raise ResourceGuardError(f"invalid Content-Length {value!r}; refusing download") from exc
    if content_length < 0:
        raise ResourceGuardError(f"invalid Content-Length {value!r}; refusing download")
    return content_length


def assess_download(content_length, partial_size, free_bytes, min_free_gb, max_file_gb):
    if min(content_length, partial_size, free_bytes) < 0:
        raise ResourceGuardError("file sizes and free space must not be negative")
    max_file_bytes = gb_to_bytes(max_file_gb)
    if content_length > max_file_bytes:
        raise ResourceGuardError(
            f"download size {format_bytes(content_length)} exceeds maximum file size "
            f"{format_bytes(max_file_bytes)}"
        )
    if partial_size > content_length:
        raise ResourceGuardError(
            f"existing partial file ({format_bytes(partial_size)}) is larger than the "
            f"remote file ({format_bytes(content_length)}); refusing to resume"
        )
    remaining_bytes = content_length - partial_size
    reserve_bytes = gb_to_bytes(min_free_gb)
    if free_bytes - remaining_bytes < reserve_bytes:
        raise ResourceGuardError(
            f"download would leave {format_bytes(free_bytes - remaining_bytes)} free; "
            f"{format_bytes(reserve_bytes)} reserve required"
        )
    return remaining_bytes


def partial_file_size(output):
    try:
        if not output.exists():
            return 0
        if not output.is_file():
            raise ResourceGuardError(f"download output is not a regular file: {output}")
        return output.stat().st_size
    except OSError as exc:
        raise ResourceGuardError(f"cannot inspect download output {output}: {exc}") from exc


def fetch_content_length(url, timeout=30):
    request = urllib.request.Request(
        url,
        method="HEAD",
        headers={"User-Agent": "agescape-resource-guard/1"},
    )
    try:
        response = urllib.request.urlopen(request, timeout=timeout)
    except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, OSError) as exc:
        raise ResourceGuardError(f"cannot retrieve download headers for {url}: {exc}") from exc
    with response:
        return parse_content_length(response.headers.get("Content-Length"))


def verify_http_range_resume(url, partial_size, existing_path=None, overlap_bytes=65536, timeout=30):
    """Verify server honors HTTP 206 Range and matches the tail of an existing partial file."""
    if partial_size <= 0:
        return True
    overlap = min(int(overlap_bytes), int(partial_size))
    start = int(partial_size) - overlap
    end = int(partial_size) - 1
    request = urllib.request.Request(
        url,
        method="GET",
        headers={
            "User-Agent": "agescape-resource-guard/1",
            "Range": f"bytes={start}-{end}",
        },
    )
    try:
        response = urllib.request.urlopen(request, timeout=timeout)
    except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, OSError) as exc:
        raise ResourceGuardError(f"range resume probe failed for {url}: {exc}") from exc
    with response:
        status = getattr(response, "status", None) or response.getcode()
        if status != 206:
            raise ResourceGuardError(
                f"server returned HTTP {status} instead of 206 Partial Content for range {start}-{end}"
            )
        remote_slice = response.read(overlap + 1)
    if len(remote_slice) != overlap:
        raise ResourceGuardError(
            f"range resume probe expected {overlap} bytes for {start}-{end}, got {len(remote_slice)}"
        )
    if existing_path is not None:
        with Path(existing_path).open("rb") as handle:
            handle.seek(start)
            local_slice = handle.read(overlap)
        if local_slice != remote_slice:
            raise ResourceGuardError(
                f"local partial tail bytes {start}-{end} do not match remote server; refusing resume"
            )
    return True


def verify_gzip_integrity(path):
    """Run gzip -t on a completed .gz download and fail closed if corrupt."""
    try:
        completed = subprocess.run(["gzip", "-t", str(path)], capture_output=True, text=True, check=False)
    except OSError as exc:
        raise ResourceGuardError(f"cannot run gzip -t on {path}: {exc}") from exc
    if completed.returncode != 0:
        detail = (completed.stderr or completed.stdout or "").strip()
        raise ResourceGuardError(f"gzip -t integrity check failed for {path}: {detail}")
    return True


def verify_min_file_size(path, min_size_bytes):
    """Verify downloaded file meets minimum byte floor."""
    actual = partial_file_size(Path(path))
    if actual < int(min_size_bytes):
        raise ResourceGuardError(
            f"downloaded file {path} is {actual} bytes, below minimum {int(min_size_bytes)} bytes"
        )
    return actual


@contextlib.contextmanager
def heavy_job_lock(lock_path=DEFAULT_LOCK_PATH):
    try:
        handle = lock_path.open("a+")
    except OSError as exc:
        raise ResourceGuardError(f"cannot open heavy-job lock {lock_path}: {exc}") from exc
    try:
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise ResourceGuardError(
                f"another heavy job holds the lock {lock_path}; refusing to start"
            ) from exc
        except OSError as exc:
            raise ResourceGuardError(f"cannot acquire heavy-job lock {lock_path}: {exc}") from exc
        yield
    finally:
        handle.close()


def run_guarded(args):
    command = build_systemd_command(args.memory_gb, normalize_command(args.command))
    with heavy_job_lock():
        check_free_space(args.path, args.min_free_gb)
        try:
            completed = subprocess.run(command, check=False)
        except OSError as exc:
            raise ResourceGuardError(f"cannot launch guarded command: {exc}") from exc
    return completed.returncode


def download_guarded(args):
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with heavy_job_lock():
        partial_size = partial_file_size(args.output)
        content_length = fetch_content_length(args.url, timeout=args.connect_timeout)
        try:
            free_bytes = shutil.disk_usage(args.output.parent).free
        except OSError as exc:
            raise ResourceGuardError(
                f"cannot check free space at {args.output.parent}: {exc}"
            ) from exc
        remaining_bytes = assess_download(
            content_length,
            partial_size,
            free_bytes,
            args.min_free_gb,
            args.max_file_gb,
        )
        if remaining_bytes == 0 and args.output.exists():
            print(f"Download already complete: {args.output}")
        else:
            if partial_size > 0 and args.verify_resume_overlap:
                verify_http_range_resume(
                    args.url,
                    partial_size,
                    existing_path=args.output,
                    timeout=args.connect_timeout,
                )
            command = build_curl_command(
                args.url,
                args.output,
                gb_to_bytes(args.max_file_gb),
                connect_timeout=args.connect_timeout,
                max_time=args.max_time,
                speed_limit=args.speed_limit,
                speed_time=args.speed_time,
            )
            try:
                completed = subprocess.run(command, check=False)
            except OSError as exc:
                raise ResourceGuardError(f"cannot launch curl: {exc}") from exc
            if completed.returncode != 0:
                return completed.returncode
        if args.min_size_bytes is not None:
            verify_min_file_size(args.output, args.min_size_bytes)
        if args.verify_gzip:
            verify_gzip_integrity(args.output)
    return 0


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="subcommand", required=True)

    run_parser = subparsers.add_parser("run", help="run a command with disk and memory guards")
    run_parser.add_argument("--path", type=Path, required=True, help="filesystem path to check")
    run_parser.add_argument("--min-free-gb", type=positive_gb, required=True)
    run_parser.add_argument("--memory-gb", "--memory-limit-gb", type=positive_gb, required=True)
    run_parser.add_argument("command", nargs=argparse.REMAINDER, help="command and arguments, after --")
    run_parser.set_defaults(handler=run_guarded)

    download_parser = subparsers.add_parser("download", help="download with disk guards and resume")
    download_parser.add_argument("--url", required=True)
    download_parser.add_argument("--output", type=Path, required=True)
    download_parser.add_argument("--min-free-gb", type=positive_gb, required=True)
    download_parser.add_argument("--max-file-gb", type=positive_gb, required=True)
    download_parser.add_argument("--connect-timeout", type=int, default=30)
    download_parser.add_argument("--max-time", type=int, default=None)
    download_parser.add_argument("--speed-limit", type=int, default=10240)
    download_parser.add_argument("--speed-time", type=int, default=120)
    download_parser.add_argument("--min-size-bytes", type=int, default=None)
    download_parser.add_argument("--verify-gzip", action="store_true")
    download_parser.add_argument("--verify-resume-overlap", action="store_true")
    download_parser.set_defaults(handler=download_guarded)
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.handler(args)
    except ResourceGuardError as exc:
        print(f"resource_guard: error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
