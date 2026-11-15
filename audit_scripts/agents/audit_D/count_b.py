import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tc_variant
buckets, stats = tc_variant.enumerate_triplane_up_to_b4(max_c_tri=int(sys.argv[1]))
print(os.environ.get("BRIDGES"), os.environ.get("ARCH_MODE"), "max_c", sys.argv[1], {k: v for k, v in stats.items() if not k.startswith("tangle") and not k.startswith("elapsed")})
