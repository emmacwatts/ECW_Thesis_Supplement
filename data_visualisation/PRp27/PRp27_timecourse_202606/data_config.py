"""Shared source-data paths for PRp27 analysis scripts.

Edit TIMECOURSE_WORKBOOK if the main Excel workbook moves again.
You can also override it temporarily with the PRP27_TIMECOURSE_WORKBOOK
environment variable, or with a script-specific variable such as GFP_SRC.
"""

from pathlib import Path
import os


TIMECOURSE_WORKBOOK = Path(
    "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/"
    "6. PosHits/prp27/prp_TimeCourse2/prp27_timecourse_Results_quantified.xlsx"
)
TIMECOURSE_WORKBOOK_NAME = "prp27_timecourse_Results_quantified.xlsx"
TIMECOURSE_WORKBOOK_ENV = "PRP27_TIMECOURSE_WORKBOOK"


def resolve_timecourse_workbook(analysis_dir=None, env_var=None):
    """Return the current main PRp27 timecourse workbook path."""
    env_vars = [env_var, TIMECOURSE_WORKBOOK_ENV] if env_var else [TIMECOURSE_WORKBOOK_ENV]
    for name in env_vars:
        env_path = os.environ.get(name)
        if env_path:
            candidate = Path(env_path).expanduser()
            if candidate.exists():
                return candidate

    if TIMECOURSE_WORKBOOK.exists():
        return TIMECOURSE_WORKBOOK

    candidate_dirs = []
    if analysis_dir is not None:
        analysis_dir = Path(analysis_dir)
        candidate_dirs.extend([analysis_dir, analysis_dir.parent])
    candidate_dirs.extend(
        [
            Path.cwd(),
            Path.cwd().parent,
            Path.home() / "Desktop" / "prp_TimeCourse2",
            Path.home() / "Desktop",
        ]
    )
    for folder in candidate_dirs:
        candidate = folder / TIMECOURSE_WORKBOOK_NAME
        if candidate.exists():
            return candidate

    raise FileNotFoundError(
        f"Could not find {TIMECOURSE_WORKBOOK_NAME}. Edit TIMECOURSE_WORKBOOK in "
        "PRp27/PRp27_timecourse_202606/data_config.py, put the workbook next to the analysis scripts, or set "
        f"{TIMECOURSE_WORKBOOK_ENV} to the workbook path."
    )
