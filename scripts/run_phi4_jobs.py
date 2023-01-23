import subprocess as proc
import time
import numpy as np

L = 16

for m2 in np.arange(1.15,1.25,0.01): 
    m2i = int(10000 * m2)
    start = time.time()
    results = proc.run(
        ['./cmake-build-release/phi4', '--lambda', "1.0", '--kappa', "1.0", '--n-term', "1000", '--n-sweeps', "1000000",
         '--save-freq',
         "10", '--seed', "34", '--Lx', f"{L}", '--Ly', f"{L}", '--M2', f"{-m2}", "--output",
         f"data/k100_m{m2i:05d}_{L:02d}x{L:02d}.bin"])
    end = time.time()
    print(f"ellapsed time {end - start:.2f}s")
