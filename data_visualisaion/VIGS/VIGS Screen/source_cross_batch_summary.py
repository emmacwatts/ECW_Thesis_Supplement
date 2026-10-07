"""Build source-derived VIGS target/batch counts from batch-folder workbooks."""

from pathlib import Path
import importlib.util

import pandas as pd


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("vigs_screen", HERE / "202607_vigs_screen.py")
VIGS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(VIGS)


def source_target(label):
    target = VIGS.canonical_line(label)
    # These are alternate labels found in the individual Batch9 workbook.
    return {"LIK1 (2)": "LIK1", "NBD": "NBD030152"}.get(target, target)


def main():
    file_rows = []
    measurement_rows = []
    batch_dirs = sorted(
        [p for p in VIGS.SOURCE_ROOT.iterdir() if p.is_dir() and p.name.lower() != "scrap"],
        key=VIGS.batch_sort_key,
    )

    for batch_dir in batch_dirs:
        files = sorted(batch_dir.rglob("*.xlsx"))
        if not files:
            file_rows.append(
                {"batch": batch_dir.name, "source_file": "", "source_sheet": "", "status": "no_xlsx", "details": ""}
            )
            continue
        for path in files:
            forced_sheet = "All" if batch_dir.name == "Batch4" and path.name == "VIGSdata_20240616.xlsx" else None
            try:
                data = VIGS.parse_excel_table(path, forced_sheet=forced_sheet)
                data["plant_line"] = data["plant_line"].map(source_target)
                if batch_dir.name == "Batch9" and "20241115" in str(path):
                    # Experimental record correction: LYK5 was not part of Batch9,
                    # despite LYK5-labelled rows being present in this workbook.
                    data = data.loc[data["plant_line"] != "LYK5"].copy()
                sheet = str(data["source_sheet"].iloc[0])
                file_rows.append(
                    {
                        "batch": batch_dir.name,
                        "source_file": str(path),
                        "source_sheet": sheet,
                        "status": "parsed",
                        "details": "Batch4 forced to All sheet" if forced_sheet else "",
                    }
                )
                for target, subset in data.groupby("plant_line"):
                    measurement_rows.append(
                        {
                            "batch": batch_dir.name,
                            "target": target,
                            "source_file": str(path),
                            "source_sheet": sheet,
                            "n_in_file": len(subset),
                        }
                    )
            except Exception as exc:
                file_rows.append(
                    {
                        "batch": batch_dir.name,
                        "source_file": str(path),
                        "source_sheet": forced_sheet or "",
                        "status": "parse_error",
                        "details": str(exc),
                    }
                )

    files = pd.DataFrame(file_rows)
    measurements = pd.DataFrame(measurement_rows)
    exclusions = VIGS.explicit_batch_exclusions()
    exclusion_lookup = {
        (row.batch, row.plant_line): row.exclusion_source for row in exclusions.itertuples()
    }
    global_exclusions = VIGS.TARGETS_EXCLUDED_FROM_SCREEN

    by_batch = (
        measurements.groupby(["batch", "target"], as_index=False)
        .agg(
            n_source=("n_in_file", "sum"),
            n_source_files=("source_file", "nunique"),
            source_files=("source_file", lambda x: "; ".join(dict.fromkeys(x))),
            source_sheets=("source_sheet", lambda x: "; ".join(dict.fromkeys(x))),
        )
    )
    reasons = []
    for row in by_batch.itertuples():
        reason = global_exclusions.get(row.target, "")
        batch_reason = exclusion_lookup.get((row.batch, row.target), "")
        reasons.append("; ".join(part for part in [reason, batch_reason] if part))
    by_batch["target_batch_exclusion_reason"] = reasons
    by_batch["target_batch_excluded"] = by_batch["target_batch_exclusion_reason"].ne("")
    by_batch["n_after_target_exclusions"] = by_batch["n_source"].where(~by_batch["target_batch_excluded"], 0)
    by_batch = by_batch.sort_values(["target", "batch"], key=lambda x: x)

    totals = (
        by_batch.groupby("target", as_index=False)
        .agg(
            n_source_across_all_batches=("n_source", "sum"),
            n_source_batches=("batch", "nunique"),
            n_after_target_exclusions=("n_after_target_exclusions", "sum"),
            n_batches_after_target_exclusions=("target_batch_excluded", lambda x: int((~x).sum())),
            n_excluded_batch_occurrences=("target_batch_excluded", "sum"),
            batches_present=("batch", lambda x: "; ".join(dict.fromkeys(x))),
        )
        .sort_values(["n_source_batches", "target"], ascending=[False, True])
    )
    breakdown = by_batch.pivot(index="target", columns="batch", values="n_source").fillna(0).astype(int)
    breakdown.insert(0, "n_source_across_all_batches", breakdown.sum(axis=1))
    breakdown.insert(1, "n_source_batches", (breakdown.drop(columns="n_source_across_all_batches") > 0).sum(axis=1))
    breakdown = breakdown.reset_index()

    batch_summary = (
        by_batch.groupby("batch", as_index=False)
        .agg(
            n_targets_in_source=("target", "nunique"),
            n_measurements_in_source=("n_source", "sum"),
            targets_in_source=("target", lambda x: "; ".join(sorted(x))),
        )
    )
    batch_summary = pd.DataFrame({"batch": [p.name for p in batch_dirs]}).merge(batch_summary, how="left", on="batch")
    batch_summary = batch_summary.merge(
        files.groupby("batch", as_index=False).agg(
            n_xlsx_records=("source_file", lambda x: int((x != "").sum())),
            file_status=("status", lambda x: "; ".join(dict.fromkeys(x))),
        ),
        how="left",
        on="batch",
    )

    output = VIGS.LOG_DIR / "vigs_screen_source_cross_batch_summary.xlsx"
    with pd.ExcelWriter(output) as writer:
        totals.to_excel(writer, sheet_name="Target totals", index=False)
        by_batch.to_excel(writer, sheet_name="Target by batch", index=False)
        breakdown.to_excel(writer, sheet_name="N matrix", index=False)
        batch_summary.to_excel(writer, sheet_name="Batch summary", index=False)
        files.to_excel(writer, sheet_name="Source files", index=False)
        measurements.to_excel(writer, sheet_name="File-level counts", index=False)
    print(output)


if __name__ == "__main__":
    main()
