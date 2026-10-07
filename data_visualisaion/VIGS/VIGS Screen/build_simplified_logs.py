"""Create three user-facing VIGS audit logs from individual batch files."""

from pathlib import Path

import pandas as pd


HERE = Path(__file__).resolve().parent
DETAIL_DIR = HERE / "202607_vigs_screen_logs"
OUTPUT_DIR = DETAIL_DIR / "current"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
SOURCE_SUMMARY = DETAIL_DIR / "vigs_screen_source_cross_batch_summary.xlsx"

QUALITY_UNIVERSAL = {
    "STP13": "cloning error",
    "GOX4": "cloning error",
    "PYL1": "cloning error",
    "LYM1": "cloning error",
    "BIK1": "universal exclusion",
    "LYK5": "universal exclusion",
    "NBD030152": "universal exclusion",
}
EXCESS_BATCH_EXCLUSIONS = [
    ("Batch1", "EAH"), ("Batch1", "EAS"), ("Batch1", "ETR1"),
    ("Batch8", "LYK3"), ("Batch1", "MPK3/6"),
    ("Batch9", "NRC234"), ("Batch1", "NRC234"), ("Batch1", "PEROX"),
    ("Batch1", "PYL1"), ("Batch5", "PYL8"),
    ("Batch1", "PYL8"), ("Batch1", "SQS"), ("Batch5", "JAZ"),
]
BIOLOGICAL_TARGETS = [
    "CDPK", "EAS", "EDS5", "ETR1", "HEVA", "JAZ",
    "PYL8", "RBOHA", "RBOHB", "SOBIR1",
]
BIOLOGICAL_UNIVERSAL = {
    "CDPK", "EAS", "EDS5", "ETR1", "HEVA", "JAZ",
    "PYL8", "RBOHA", "RBOHB", "SOBIR1",
}
FAILED_SCREEN_ROUNDS = {
    "Batch12": ["ICS", "EAS", "ETR1"],
    "Batch13": ["PYL1", "PYL8"],
}


def joined(values):
    return "; ".join(sorted({str(value) for value in values if pd.notna(value) and str(value)}))


def main():
    by_batch = pd.read_excel(SOURCE_SUMMARY, sheet_name="Target by batch")
    batch_summary = pd.read_excel(SOURCE_SUMMARY, sheet_name="Batch summary")
    source_files = pd.read_excel(SOURCE_SUMMARY, sheet_name="Source files")
    manifest = pd.read_excel(DETAIL_DIR / "vigs_screen_batch_manifest.xlsx")
    normalized = pd.read_excel(DETAIL_DIR / "vigs_screen_normalized_long.xlsx")
    batch_summary["nonquantitative_screen_outcomes"] = batch_summary["batch"].map(
        {batch: "; ".join(targets) for batch, targets in FAILED_SCREEN_ROUNDS.items()}
    )
    batch_summary["outcome_note"] = batch_summary["batch"].map(
        {
            batch: "all listed candidates failed; count as screen attempts, no quantitative GFP data"
            for batch in FAILED_SCREEN_ROUNDS
        }
    )
    included_by_batch = (
        normalized.loc[normalized["plant_line"] != "GUS"]
        .groupby("batch")["plant_line"]
        .apply(lambda x: "; ".join(sorted(set(x))))
        .to_dict()
    )
    biological_excluded_by_batch = {}
    for batch, subset in manifest.loc[manifest["manifest_status"] == "excluded"].groupby("batch"):
        targets = {row.plant_line for row in subset.itertuples() if row.plant_line not in QUALITY_UNIVERSAL}
        targets |= {target for target in FAILED_SCREEN_ROUNDS.get(batch, []) if target not in QUALITY_UNIVERSAL}
        biological_excluded_by_batch[batch] = "; ".join(sorted(targets))
    for batch, targets in FAILED_SCREEN_ROUNDS.items():
        existing = set(filter(None, biological_excluded_by_batch.get(batch, "").split("; ")))
        nonquality_targets = {target for target in targets if target not in QUALITY_UNIVERSAL}
        biological_excluded_by_batch[batch] = "; ".join(sorted(existing | nonquality_targets))
    quality_by_batch = {}
    for target in QUALITY_UNIVERSAL:
        source_batches = set(by_batch.loc[by_batch["target"] == target, "batch"].astype(str))
        recorded_batches = set(manifest.loc[manifest["plant_line"] == target, "batch"].astype(str))
        for batch in source_batches | recorded_batches:
            quality_by_batch.setdefault(batch, set()).add(target)
    for batch, target in EXCESS_BATCH_EXCLUSIONS:
        quality_by_batch.setdefault(batch, set()).add(f"{target} (excess batch)")
    quality_by_batch.setdefault("Batch3", set()).update(
        {"PAL (second-date set)", "STP13 (second-date set)"}
    )
    quality_by_batch.setdefault("Batch9", set()).update(
        {
            "EAH (20241115)", "LYK3 (20241115)",
            "NRC234 (20241115)", "SQS (20241115)", "LIK1 (20241217 partial set)",
        }
    )
    batch_summary["included_targets_results"] = batch_summary["batch"].map(included_by_batch)
    batch_summary["excluded_targets_results"] = batch_summary["batch"].map(biological_excluded_by_batch)
    batch_summary["quality_excluded_targets"] = batch_summary["batch"].map(
        {batch: "; ".join(sorted(targets)) for batch, targets in quality_by_batch.items()}
    )

    with pd.ExcelWriter(OUTPUT_DIR / "01_batch_metadata.xlsx") as writer:
        batch_summary.to_excel(writer, sheet_name="Batch summary", index=False)
        source_files.to_excel(writer, sheet_name="Individual source files", index=False)
        by_batch[["batch", "target", "n_source", "source_files", "source_sheets"]].to_excel(
            writer, sheet_name="Target n by batch", index=False
        )

    quality_rows = []
    for target, reason in QUALITY_UNIVERSAL.items():
        source = by_batch.loc[by_batch["target"] == target]
        quality_rows.append(
            {
                "target": target,
                "scope": "universal",
                "reason": reason,
                "affected_batches_in_source": joined(source["batch"]),
                "n_removed_from_source": int(source["n_source"].sum()),
                "note": "",
            }
        )
    batch3_second_date = by_batch.loc[
        (by_batch["batch"] == "Batch3") & by_batch["target"].isin(["PAL", "STP13"])
    ]
    for target in ["PAL", "STP13"]:
        source = batch3_second_date.loc[batch3_second_date["target"] == target]
        quality_rows.append(
            {
                "target": target,
                "scope": "Batch3 second-date set",
                "reason": "no same-day GUS control",
                "affected_batches_in_source": "Batch3",
                "n_removed_from_source": int(source["n_source"].sum()),
                "note": "STP13 is also universally excluded for cloning error" if target == "STP13" else "",
            }
        )
    quality_rows.extend(
        [
            {
                "target": "EAH; LYK3; NRC234; SQS",
                "scope": "Batch9 20241115 set",
                "reason": "no same-day GUS control",
                "affected_batches_in_source": "Batch9",
                "n_removed_from_source": 40,
                "note": "entire 20241115 workbook excluded from quantitative normalization",
            },
            {
                "target": "LIK1",
                "scope": "Batch9 20241217 set",
                "reason": "no same-day GUS control",
                "affected_batches_in_source": "Batch9",
                "n_removed_from_source": 16,
                "note": "20241204 LIK1 n=6 retained with same-day GUS",
            },
        ]
    )
    for batch, target in EXCESS_BATCH_EXCLUSIONS:
        source = by_batch.loc[(by_batch["batch"] == batch) & (by_batch["target"] == target)]
        quality_rows.append(
            {
                "target": target,
                "scope": batch,
                "reason": "excess batch; only three highest-numbered batch attempts retained",
                "affected_batches_in_source": batch,
                "n_removed_from_source": int(source["n_source"].sum()),
                "note": "ordered by numeric batch sequence",
            }
        )
    pd.DataFrame(quality_rows).to_excel(OUTPUT_DIR / "02_quality_exclusions.xlsx", index=False)

    biological_rows = []
    for target in BIOLOGICAL_TARGETS:
        source = by_batch.loc[by_batch["target"] == target]
        recorded = manifest.loc[manifest["plant_line"] == target]
        excluded = recorded.loc[recorded["manifest_status"] == "excluded"]
        source_batches = set(source["batch"].dropna().astype(str))
        excluded_batches = set(excluded["batch"].dropna().astype(str))
        retained = source.loc[~source["batch"].astype(str).isin(excluded_batches)]
        biological_rows.append(
            {
                "target": target,
                "universal_exclusion_by_at_least_two_failed_batches": target in BIOLOGICAL_UNIVERSAL,
                "source_batches": joined(source_batches),
                "n_source_total": int(source["n_source"].sum()),
                "batches_recorded_as_excluded": joined(excluded_batches),
                "excluded_batches_with_source_data": joined(source_batches & excluded_batches),
                "n_removed_where_source_exists": int(
                    source.loc[source["batch"].astype(str).isin(excluded_batches), "n_source"].sum()
                ),
                "retained_source_batches": joined(retained["batch"]),
                "n_retained": int(retained["n_source"].sum()),
            }
        )
    failed_rows = [
        {
            "batch": batch,
            "target": target,
            "screen_outcome": "failed",
            "counts_as_batch_attempt": True,
            "quantitative_gfp_available": False,
            "exclusion_type": "screen result (not quality-based)",
        }
        for batch, targets in FAILED_SCREEN_ROUNDS.items()
        for target in targets
        if target not in QUALITY_UNIVERSAL
    ]
    with pd.ExcelWriter(OUTPUT_DIR / "03_biological_exclusions.xlsx") as writer:
        pd.DataFrame(biological_rows).to_excel(writer, sheet_name="Phenotype exclusions", index=False)
        pd.DataFrame(failed_rows).to_excel(writer, sheet_name="Failed screen rounds", index=False)

    failed = pd.DataFrame(failed_rows)
    candidates = sorted(
        (set(by_batch["target"]) | set(manifest["plant_line"]) | set(failed["target"])) - {"GUS"}
    )
    candidate_rows = []
    for target in candidates:
        source = by_batch.loc[by_batch["target"] == target]
        recorded = manifest.loc[manifest["plant_line"] == target]
        failed_target = failed.loc[failed["target"] == target]
        retained = normalized.loc[normalized["plant_line"] == target]
        attempted_batches = (
            set(source["batch"].dropna().astype(str))
            | set(recorded["batch"].dropna().astype(str))
            | set(failed_target["batch"].dropna().astype(str))
        )
        biological_excluded = (
            set(recorded.loc[recorded["manifest_status"] == "excluded", "batch"].astype(str))
            | set(failed_target["batch"].astype(str))
        )
        if target in QUALITY_UNIVERSAL:
            biological_excluded = set()
        quality_excluded = set()
        fully_quality_excluded_batches = set()
        if target in QUALITY_UNIVERSAL:
            quality_excluded |= attempted_batches
            fully_quality_excluded_batches |= attempted_batches
        if target in {"PAL", "STP13"}:
            quality_excluded.add("Batch3")
            fully_quality_excluded_batches.add("Batch3")
        if target in {"EAH", "LYK3", "NRC234"}:
            quality_excluded.add("Batch9 (20241115)")
            fully_quality_excluded_batches.add("Batch9")
        if target == "SQS":
            quality_excluded.add("Batch9 (20241115)")
        if target == "LIK1":
            quality_excluded.add("Batch9 (20241217 partial set)")
        excess_batches = {batch for batch, excess_target in EXCESS_BATCH_EXCLUSIONS if excess_target == target}
        quality_excluded |= excess_batches
        fully_quality_excluded_batches |= excess_batches
        retained_batches = set(retained["batch"].dropna().astype(str))
        source_batches = set(source["batch"].dropna().astype(str))
        listed_included = set(
            recorded.loc[recorded["manifest_status"] == "included", "batch"].dropna().astype(str)
        )
        listed_included_without_source = listed_included - source_batches
        n_by_batch = "; ".join(
            f"{row.batch}: n={int(row.n_source)}" for row in source.sort_values("batch").itertuples()
        )
        retained_n_by_batch = "; ".join(
            f"{batch}: n={len(group)}" for batch, group in retained.groupby("batch", sort=True)
        )
        valid_batches = attempted_batches - fully_quality_excluded_batches
        candidate_rows.append(
            {
                "target": target,
                "batches_recorded_or_present": joined(attempted_batches),
                "n_batches_recorded_or_present": len(attempted_batches),
                "n_valid_batches": len(valid_batches),
                "source_n_by_batch": n_by_batch,
                "source_n_total": int(source["n_source"].sum()),
                "included_batches": joined(retained_batches),
                "times_included": len(retained_batches),
                "listed_included_but_no_source_batches": joined(listed_included_without_source),
                "times_listed_included_but_no_source": len(listed_included_without_source),
                "retained_n_by_batch": retained_n_by_batch,
                "retained_n_total": len(retained),
                "biological_or_failed_screen_exclusions": joined(biological_excluded),
                "times_biologically_excluded_or_failed": len(biological_excluded),
                "quality_exclusions": joined(quality_excluded),
                "times_quality_excluded": len(quality_excluded),
            }
        )
    pd.DataFrame(candidate_rows).to_excel(OUTPUT_DIR / "04_candidate_batch_summary.xlsx", index=False)
    print(OUTPUT_DIR)


if __name__ == "__main__":
    main()
