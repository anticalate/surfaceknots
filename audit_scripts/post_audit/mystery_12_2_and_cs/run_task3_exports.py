"""Regenerate synchronized family_catalog.{json,csv}, hub_realizations.json,
complexity_comparison.{json,csv}, and surfaceknots_atlas.html with strictly <= 8 workers.
"""
from __future__ import annotations

import os
import sys
import time
from pathlib import Path

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

REPO = Path("/afs/csail.mit.edu/u/e/ejli/.gemini/antigravity-cli/worktrees/1dadd8f0-f304-455a-b682-76a831212ce6/subagent-4D-Topology---Mystery-Knot-Analyst-self-69be3150")
sys.path.insert(0, str(REPO / "code" / "src"))

from surfaceknots.atlas import export_interactive_atlas
from surfaceknots.complexity import export_complexity_comparison_datasets
from surfaceknots.families import export_family_datasets


def main() -> None:
    data_dir = REPO / "code" / "data"
    t0 = time.perf_counter()
    print("1. Running export_family_datasets (num_workers=6, <= 8 cores)...")
    fam_res = export_family_datasets(data_dir, num_workers=6)
    print(f"   Finished in {time.perf_counter() - t0:.2f}s!")
    print(
        f"   hub_max_ch={fam_res['hub_max_ch']}, total_instances={fam_res['total_family_instances']}, "
        f"hub_matched_instances={fam_res['hub_matched_instances']}, "
        f"distinct_hub_realized={fam_res['distinct_hub_entries_realized']}, "
        f"distinct_new_spoke_classes={fam_res['distinct_ch_ge_12_invariant_classes']}"
    )
    print(f"   exact_ch12={fam_res['exact_ch_12_certified_instances']}")
    print(f"   exact_ch13={fam_res['exact_ch_13_certified_instances']}")

    t1 = time.perf_counter()
    print("2. Running export_complexity_comparison_datasets...")
    comp_res = export_complexity_comparison_datasets(data_dir)
    print(f"   Finished in {time.perf_counter() - t1:.2f}s! total_distinct_classes_compared={comp_res['total_distinct_classes_compared']}")

    t2 = time.perf_counter()
    print("3. Running export_interactive_atlas...")
    atlas_out = REPO / "expo" / "surfaceknots_atlas.html"
    atlas_counts = export_interactive_atlas(data_dir, atlas_out)
    print(f"   Finished in {time.perf_counter() - t2:.2f}s! atlas_counts={atlas_counts}")


if __name__ == "__main__":
    main()
