#!/usr/bin/env python3
from pathlib import Path
import pickle
import sys
import time

REPO_ROOT = Path("/afs/csail.mit.edu/u/e/ejli/.gemini/antigravity-cli/worktrees/1dadd8f0-f304-455a-b682-76a831212ce6/subagent-C99---Unknotter-Speedup-Engineer-self-1528875b")
SRC_DIR = REPO_ROOT / "code" / "src"
SCRATCH_DIR = Path("/afs/csail.mit.edu/u/e/ejli/.gemini/antigravity-cli/brain/ae1d4843-4647-4011-8a0b-eb0a15e61a22/scratch")
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from surfaceknots.complexity import analyze_connected_sum_collapses, build_complexity_comparison
from surfaceknots.families import build_family_catalog


def main() -> None:
    data_dir = REPO_ROOT / "code" / "data"

    t0 = time.perf_counter()
    col = analyze_connected_sum_collapses()
    dt_col = time.perf_counter() - t0
    print(f"analyze_connected_sum_collapses: {dt_col:.3f}s ({len(col)} instances)")

    t0 = time.perf_counter()
    cmp_res = build_complexity_comparison(data_dir)
    dt_cmp = time.perf_counter() - t0
    print(f"build_complexity_comparison:     {dt_cmp:.3f}s ({cmp_res['total_distinct_classes_compared']} classes)")

    t0 = time.perf_counter()
    fam_res = build_family_catalog(data_dir)
    dt_fam = time.perf_counter() - t0
    print(f"build_family_catalog:            {dt_fam:.3f}s ({fam_res['total_family_instances']} instances)")

    with open(SCRATCH_DIR / "oracle_datasets.pkl", "wb") as f:
        pickle.dump(
            {"col": col, "cmp_res": cmp_res, "fam_res": fam_res, "timings": {"col": dt_col, "cmp": dt_cmp, "fam": dt_fam}},
            f,
        )
    print("Saved oracle_datasets.pkl")


if __name__ == "__main__":
    main()
