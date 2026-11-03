"""Portable paths shared by analysis entry points."""
from pathlib import Path
import os
REPO_ROOT = Path(os.environ.get("IL10_PROJECT_ROOT", Path(__file__).resolve().parents[1])).resolve()
DATA_DIR = Path(os.environ.get("IL10_DATA_DIR", REPO_ROOT / "data")).expanduser().resolve()
RESULTS_DIR = Path(os.environ.get("IL10_RESULTS_DIR", REPO_ROOT / "outputs")).expanduser().resolve()
ANNOTATION_DIR = Path(os.environ.get("IL10_ANNOTATION_DIR", DATA_DIR / "annotation")).expanduser().resolve()
def data_path(*parts): return DATA_DIR.joinpath(*parts)
def result_path(*parts): return RESULTS_DIR.joinpath(*parts)
def annotation_path(*parts): return ANNOTATION_DIR.joinpath(*parts)
