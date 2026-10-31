from pathlib import Path

import pandas as pd


def load_program_tiers(repo_root: Path, target: str = "human") -> dict[str, dict[str, list[str]]]:
    repo_root = Path(repo_root)
    hierarchy = pd.read_csv(repo_root / "results" / "m5-program-hierarchy.csv")
    if target == "mouse":
        mapper = {}
    elif target == "human":
        orthologs = pd.read_csv(repo_root / "data" / "reference" / "orthologs-mouse-human.csv")
        mapper = dict(zip(orthologs["mouse_gene_symbol"], orthologs["human_gene_symbol"]))
    else:
        raise ValueError(f"Unsupported ortholog target: {target}")

    programs = {}
    for tier in ("P0_Raw257", "P1_Clean"):
        tier_rows = hierarchy[hierarchy["tier"] == tier]
        if tier_rows.empty:
            raise ValueError(f"Missing signature tier: {tier}")
        programs[tier] = {}
        for direction, key in (("UP", "up"), ("DOWN", "down")):
            mouse_genes = tier_rows.loc[tier_rows["direction"] == direction, "gene"].tolist()
            if target == "mouse":
                genes = mouse_genes
            else:
                genes = [mapper[g] for g in mouse_genes if isinstance(mapper.get(g), str)]
            programs[tier][key] = sorted(set(genes))
    return programs
