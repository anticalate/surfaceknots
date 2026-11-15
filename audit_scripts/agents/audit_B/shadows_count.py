import time, ctypes, sys
from surfaceknots.fast_enum import _get_lib, generate_prime_shadows_c_list
from surfaceknots.quadrangulations import generate_quadrangulations, quadrangulation_to_shadow
lib = _get_lib()
maxn = int(sys.argv[1])
for n in range(2, maxn+1):
    t=time.time()
    c = int(lib.generate_prime_shadows_c(n, None, 0))
    print(n, "C count", c, "%.1fs"%(time.time()-t), flush=True)
