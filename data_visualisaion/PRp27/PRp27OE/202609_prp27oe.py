"""Replot the PRp27 overexpression dataset in the PRp27-1 mutant background."""

from pathlib import Path
from zipfile import ZipFile
import csv
import io
import os

MPLCONFIGDIR = Path.cwd() / ".matplotlib-cache"
MPLCONFIGDIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIGDIR))
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


HERE = Path(__file__).resolve().parent
SOURCE = Path(
    "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/6. PosHits/"
    "prp27/prpOE/20250418_prp27 exp/20250426_Prism_PBKM.prism"
)
FIG_DIR = HERE / "figures"
LOG_DIR = HERE / "logs"
FIG_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)

TABLE = "data/tables/131E1965-E5D1-48E9-BC06-C701ED66C520/data.csv"
GROUPS = ["EV", "PRp27 OE", "PRp27 OE E123Q"]
COLORS = ["#4E79A7", "#F28E2B", "#59A14F"]
MARKERS = ["o", "^", "s"]


def load_prism_table():
    with ZipFile(SOURCE) as archive:
        rows = list(csv.reader(io.StringIO(archive.read(TABLE).decode("utf-8"))))
    values = np.asarray(rows, dtype=float)
    return pd.DataFrame(
        {
            "background": "prp27-1",
            "group": np.repeat(GROUPS, values.shape[0]),
            "biological_replicate": np.tile(np.arange(1, values.shape[0] + 1), values.shape[1]),
            "gfp_signal": values.T.reshape(-1),
        }
    )


data = load_prism_table()
data.to_csv(LOG_DIR / "prp27oe_prp27_1_background_values.csv", index=False)

# Two planned comparisons against EV. Replicates are independent and the source
# variances are unequal, so use two-sided Welch tests with Holm correction.
ev = data.loc[data["group"] == "EV", "gfp_signal"].to_numpy()
comparison_rows = []
for group in GROUPS[1:]:
    treatment = data.loc[data["group"] == group, "gfp_signal"].to_numpy()
    result = stats.ttest_ind(ev, treatment, equal_var=False, alternative="two-sided")
    comparison_rows.append(
        {
            "comparison": f"EV vs {group}",
            "test": "two-sided Welch t-test",
            "t": result.statistic,
            "df": result.df,
            "raw_p": result.pvalue,
        }
    )
adjusted = multipletests([row["raw_p"] for row in comparison_rows], method="holm")[1]
for row, adjusted_p in zip(comparison_rows, adjusted):
    row["holm_adjusted_p"] = adjusted_p
    row["summary"] = "*" if adjusted_p < 0.05 else "ns"
comparisons = pd.DataFrame(comparison_rows)
comparisons.to_csv(LOG_DIR / "prp27oe_welch_holm_comparisons.csv", index=False)

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 14,
        "axes.labelsize": 14,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "axes.edgecolor": "#4A4A4A",
        "axes.labelcolor": "#4A4A4A",
        "xtick.color": "#4A4A4A",
        "ytick.color": "#4A4A4A",
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
    }
)

fig, ax = plt.subplots(figsize=(7.4, 4.7))
ax.set_facecolor("#E9EFF6")
ax.grid(axis="y", color="white", linewidth=1.2, zorder=0)
x = np.arange(len(GROUPS))
tops = []
for i, group in enumerate(GROUPS):
    values = data.loc[data["group"] == group, "gfp_signal"].to_numpy() / 1000
    mean = values.mean()
    sem = values.std(ddof=1) / np.sqrt(len(values))
    tops.append(max(values.max(), mean + sem))
    ax.bar(i, mean, width=0.58, color=COLORS[i], alpha=0.60, edgecolor="#303030", linewidth=1.2, zorder=2)
    ax.errorbar(i, mean, yerr=sem, color="#303030", capsize=4, linewidth=1.2, zorder=3)
    jitter = np.linspace(-0.10, 0.10, len(values))
    ax.scatter(i + jitter, values, s=44, marker=MARKERS[i], color=COLORS[i],
               alpha=0.90, edgecolor="#000000", linewidth=0.45, zorder=4)

def comparison_line(left, right, y, label):
    ax.plot([left, right], [y, y], color="#303030", linewidth=1.2, clip_on=False)
    ax.text((left + right) / 2, y + 0.18, label, ha="center", va="bottom", color="#303030", fontsize=13)

comparison_line(0, 1, max(tops) + 0.7, comparisons.iloc[0]["summary"])
comparison_line(0, 2, max(tops) + 2.0, comparisons.iloc[1]["summary"])
ax.set_ylim(-0.5, max(tops) + 4.0)
ax.set_xticks(x)
ax.set_xticklabels(["EV", "PRp27 OE", "PRp27 OE\nE123Q"])
ax.set_ylabel("GFP signal intensity (x1000)")
ax.spines[["top", "right"]].set_visible(False)
fig.subplots_adjust(left=0.14, right=0.97, top=0.95, bottom=0.20)

for extension in ("png", "svg", "pdf"):
    fig.savefig(FIG_DIR / f"202609_prp27oe_prp27_1_background.{extension}", dpi=300, bbox_inches="tight")
plt.close(fig)

(LOG_DIR / "methods_note.txt").write_text(
    "Values were read directly from the prp27- Prism data table in 20250426_Prism_PBKM.prism. "
    "The separate WT-background table was excluded as requested. Bars show mean +/- SEM and "
    "points show the four independent biological observations. Display values are divided by 1000. "
    "The two prespecified control comparisons (EV vs PRp27 OE and EV vs PRp27 OE E123Q) use "
    "independent-sample, two-sided Welch t-tests because variances are unequal. P-values are "
    "Holm-adjusted across the two comparisons. See prp27oe_welch_holm_comparisons.csv.\n"
)

print(f"Saved figure to {FIG_DIR}")
