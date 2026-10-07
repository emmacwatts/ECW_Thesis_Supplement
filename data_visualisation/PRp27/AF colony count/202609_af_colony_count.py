"""Replot the 20260420 AF colony-count experiment in thesis house style."""

from pathlib import Path
from zipfile import ZipFile
import csv
import io
from itertools import combinations, product
import os

MPLCONFIGDIR = Path.cwd() / ".matplotlib-cache"
MPLCONFIGDIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIGDIR))
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd


HERE = Path(__file__).resolve().parent
SOURCE = Path(
    "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/6. PosHits/"
    "prp27/AF_GrowthAssay/20260420/20260420_Agro_PRp27_ColonyCount.prism"
)
RAW_TABLE = "data/tables/DFD54C58-E0CD-4585-B769-384F562EDD01/data.csv"
FIG_DIR = HERE / "figures"
LOG_DIR = HERE / "logs"
FIG_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)

GROUPS = ["Bacteria only", "Agro + EV AF", "Agro + PRp27 AF"]
COLORS = ["#4E79A7", "#F28E2B", "#59A14F"]
MARKERS = ["o", "^", "s"]
TEXT_COLOR = "#4A4A4A"
ALPHA = 0.05


def load_prism_data():
    with ZipFile(SOURCE) as archive:
        rows = list(csv.reader(io.StringIO(archive.read(RAW_TABLE).decode("utf-8"))))
    rows = [row for row in rows if row and row[0].strip()]
    records = []
    for row in rows:
        time_h = float(row[0])
        for group_index, group in enumerate(GROUPS):
            for replicate in range(3):
                value = float(row[1 + group_index * 3 + replicate])
                records.append(
                    {
                        "time_h": time_h,
                        "group": group,
                        "replicate": replicate + 1,
                        "cfu_per_ml": value,
                        "log10_cfu_per_ml": np.log10(value),
                    }
                )
    return pd.DataFrame(records)


def style_axis(ax):
    ax.set_facecolor("#E9EFF6")
    ax.grid(axis="y", color="white", linewidth=1.2, zorder=0)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(colors=TEXT_COLOR, width=1.0, length=4.5)
    ax.set_axisbelow(True)


def run_timepoint_statistics(frame):
    """Compare all three treatments separately at each time on the log10 scale."""
    anova_rows = []
    tukey_rows = []
    for time_h in sorted(frame["time_h"].unique()):
        at_time = frame.loc[frame["time_h"].eq(time_h)]
        samples = [
            at_time.loc[at_time["group"].eq(group), "log10_cfu_per_ml"].to_numpy()
            for group in GROUPS
        ]
        result = stats.f_oneway(*samples)
        anova_rows.append(
            {
                "time_h": time_h,
                "response": "log10_cfu_per_ml",
                "f_statistic": result.statistic,
                "p_value": result.pvalue,
            }
        )
        tukey = pairwise_tukeyhsd(
            at_time["log10_cfu_per_ml"], at_time["group"], alpha=ALPHA
        )
        for row in tukey.summary().data[1:]:
            tukey_rows.append(
                {
                    "time_h": time_h,
                    "group_1": row[0],
                    "group_2": row[1],
                    "mean_difference": row[2],
                    "p_adjusted": row[3],
                    "ci_lower": row[4],
                    "ci_upper": row[5],
                    "significant": bool(row[6]),
                }
            )
    return pd.DataFrame(anova_rows), pd.DataFrame(tukey_rows)


def compact_letters(tukey_rows, means):
    """Return a minimal compact-letter display for three Tukey comparisons."""
    significant = {
        frozenset((row.group_1, row.group_2)): bool(row.significant)
        for row in tukey_rows.itertuples()
    }
    alphabet = "abc"
    masks = range(1, 2 ** len(alphabet))
    candidates = []
    for mask_values in product(masks, repeat=len(GROUPS)):
        assignment = {
            group: {alphabet[i] for i in range(len(alphabet)) if mask & (1 << i)}
            for group, mask in zip(GROUPS, mask_values)
        }
        valid = all(
            bool(assignment[left] & assignment[right])
            != significant.get(frozenset((left, right)), False)
            for left, right in combinations(GROUPS, 2)
        )
        if valid:
            used = set().union(*assignment.values())
            candidates.append((len(used), sum(map(len, assignment.values())), mask_values, assignment))

    if not candidates:
        return {group: alphabet[index] for index, group in enumerate(GROUPS)}

    assignment = min(candidates, key=lambda item: item[:3])[3]
    # Canonicalise so the highest mean receives "a", followed by b, c as needed.
    remap = {}
    for group in means.sort_values(ascending=False).index:
        for letter in sorted(assignment[group]):
            if letter not in remap:
                remap[letter] = alphabet[len(remap)]
    return {group: "".join(sorted(remap[letter] for letter in assignment[group])) for group in GROUPS}


def add_timecourse_letters(ax, frame, tukey_table):
    """Stack treatment-coloured compact letters over every timepoint."""
    y_min, y_max = ax.get_ylim()
    y_range = y_max - y_min
    gap = 0.035 * y_range
    step = 0.055 * y_range
    highest_label = y_max
    for time_h in sorted(frame["time_h"].unique()):
        at_time = frame.loc[frame["time_h"].eq(time_h)]
        means = at_time.groupby("group")["log10_cfu_per_ml"].mean().reindex(GROUPS)
        rows = tukey_table.loc[tukey_table["time_h"].eq(time_h)]
        letters = compact_letters(rows, means)
        base = at_time["log10_cfu_per_ml"].max() + gap
        ranked = means.sort_values(ascending=True).index
        for level, group in enumerate(ranked):
            y = base + level * step
            group_index = GROUPS.index(group)
            ax.text(
                time_h, y, letters[group], color=COLORS[group_index], fontsize=12,
                fontweight="bold", ha="center", va="bottom", zorder=7,
            )
            highest_label = max(highest_label, y + step)
    ax.set_ylim(y_min, highest_label)


data = load_prism_data()
data.to_csv(LOG_DIR / "af_colony_count_values.csv", index=False)

zero = data[data["time_h"] == 0].copy()
anova_table, tukey_table = run_timepoint_statistics(data)
anova_table.to_csv(LOG_DIR / "af_colony_count_timepoint_anova.csv", index=False)
tukey_table.to_csv(LOG_DIR / "af_colony_count_timepoint_tukey.csv", index=False)

# Retain the original 0-h filenames for downstream references, now using the
# same log10-scale analysis as the time-course panel.
zero_anova = anova_table.loc[anova_table["time_h"].eq(0)].iloc[0]
tukey_table.loc[tukey_table["time_h"].eq(0)].to_csv(
    LOG_DIR / "af_colony_count_0h_tukey.csv", index=False
)
(LOG_DIR / "af_colony_count_0h_anova.txt").write_text(
    "Ordinary one-way ANOVA of log10 CFU/mL at 0 h: "
    f"F={zero_anova['f_statistic']:.6g}, p={zero_anova['p_value']:.6g}.\n"
)

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 14,
        "axes.labelsize": 14,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "legend.fontsize": 13,
        "text.color": TEXT_COLOR,
        "axes.labelcolor": TEXT_COLOR,
        "axes.edgecolor": TEXT_COLOR,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
    }
)

fig, (ax_a, ax_b) = plt.subplots(
    1,
    2,
    figsize=(11.6, 5.2),
    gridspec_kw={"width_ratios": [1.35, 1.0]},
)

# Original time-course structure, displayed on the analysed log10 scale.
for index, group in enumerate(GROUPS):
    subset = data[data["group"] == group]
    summary = subset.groupby("time_h")["log10_cfu_per_ml"].agg(["mean", "sem"]).reindex([0, 2, 8])
    times = summary.index.to_numpy(dtype=float)
    means = summary["mean"].to_numpy()
    sems = summary["sem"].to_numpy()
    ax_a.fill_between(times, means - sems, means + sems, color=COLORS[index], alpha=0.18, linewidth=0, zorder=1)
    ax_a.plot(times, means, color=COLORS[index], marker=MARKERS[index], markeredgecolor="#303030", linewidth=2.0, markersize=7, label=group, zorder=3)
    for time_h in times:
        values = subset.loc[subset["time_h"] == time_h, "log10_cfu_per_ml"].to_numpy()
        jitter = np.linspace(-0.12, 0.12, len(values))
        ax_a.scatter(
            time_h + jitter,
            values,
            s=36,
            marker=MARKERS[index],
            color=COLORS[index],
            alpha=0.90,
            edgecolor="#000000",
            linewidth=0.45,
            zorder=4,
        )
ax_a.set_xticks([0, 2, 8])
ax_a.set_xlabel("Time (h)")
ax_a.set_ylabel(r"Colony count (log$_{10}$ CFU mL$^{-1}$)")
style_axis(ax_a)
add_timecourse_letters(ax_a, data, tukey_table)

# 0-hour raw values, scaled for readable tick labels.
x_positions = np.arange(len(GROUPS))
for index, group in enumerate(GROUPS):
    values = zero.loc[zero["group"] == group, "cfu_per_ml"].to_numpy() / 1e9
    mean = values.mean()
    sem = stats.sem(values)
    ax_b.bar(index, mean, width=0.58, color=COLORS[index], alpha=0.60, edgecolor="#303030", linewidth=1.2, zorder=2)
    ax_b.errorbar(index, mean, yerr=sem, color="#303030", capsize=4, linewidth=1.2, zorder=3)
    ax_b.scatter(
        index + np.linspace(-0.10, 0.10, len(values)),
        values,
        s=42,
        marker=MARKERS[index],
        color=COLORS[index],
        alpha=0.90,
        edgecolor="#000000",
        linewidth=0.45,
        zorder=4,
    )
ax_b.set_xticks(x_positions)
ax_b.set_xticklabels(["Bacteria\nonly", "Agro +\nEV AF", "Agro +\nPRp27 AF"])
ax_b.set_ylabel(r"Colony count (×10$^9$ CFU mL$^{-1}$)")
style_axis(ax_b)

zero_means = zero.groupby("group")["log10_cfu_per_ml"].mean().reindex(GROUPS)
zero_letters = compact_letters(tukey_table.loc[tukey_table["time_h"].eq(0)], zero_means)
right_top = ax_b.get_ylim()[1]
letter_gap = 0.035 * right_top
for index, group in enumerate(GROUPS):
    values = zero.loc[zero["group"].eq(group), "cfu_per_ml"].to_numpy() / 1e9
    y = max(values.max(), values.mean() + stats.sem(values)) + letter_gap
    ax_b.text(
        index, y, zero_letters[group], color=COLORS[index], fontsize=12,
        fontweight="bold", ha="center", va="bottom", zorder=7,
    )
ax_b.set_ylim(0, max(ax_b.get_ylim()[1], max(
    zero.loc[zero["group"].eq(group), "cfu_per_ml"].max() / 1e9 + 3 * letter_gap
    for group in GROUPS
)))

fig.legend(
    handles=ax_a.get_legend_handles_labels()[0],
    labels=GROUPS,
    loc="upper center",
    bbox_to_anchor=(0.53, 0.995),
    ncol=3,
    frameon=False,
)
fig.subplots_adjust(left=0.09, right=0.98, top=0.84, bottom=0.18, wspace=0.34)

for extension in ("png", "svg", "pdf"):
    fig.savefig(FIG_DIR / f"202609_af_colony_count_timecourse_0h.{extension}", dpi=300, bbox_inches="tight")
plt.close(fig)
. Compact letters encode the "
    "Tukey results: groups sharing a letter are not significantly different; groups with no letter in common are "
    "significantly different. The right panel displays the same 0-h compact-letter result as the left panel.\n"
)

print(f"Loaded {len(data)} observations from {SOURCE.name}")
print(f"0-h ANOVA: F={zero_anova['f_statistic']:.4g}, p={zero_anova['p_value']:.4g}")
print(f"Saved figures to: {FIG_DIR}")
