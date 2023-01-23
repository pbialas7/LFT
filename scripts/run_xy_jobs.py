import subprocess as proc
import time
import numpy as np

L = 64

for beta in np.arange(2.0/3.0,3.0,0.2):
    b2i = int(1000 * beta)
    start = time.time()
    results = proc.run(
        ['./cmake-build-release/xy', '--n-term', "1000", '--n-sweeps', "40000",
         '--seed', "34", '--Lx', f"{L}", '--Ly', f"{L}", '--beta', f"{beta}",
         "-m", "10",
         "--meas-file-name", f"data/am_b{b2i:05d}_{L:03d}x{L:03d}.txt"])
    end = time.time()
    print(f"ellapsed time {end - start:.2f}s")
