import subprocess as proc
import time
import numpy as np

L = 8
m2 = -4

n_sweeps = 4000000

for lamda in (27,):
    lamda = float(lamda)
    start = time.time()
    results = proc.run(
        ['./cmake-build-release/phi4', '--lambda', f"{lamda}", '--kappa', "1.0", '--n-term', "1000", '--n-sweeps',
         f"{n_sweeps}",
         '--save-freq',
         "10", '--seed', "34", '--Lx', f"{L}", '--Ly', f"{L}", '--M2', f"{m2}", "--output",
         f"phi4_data/cfgs_{L:02d}x{L:02d}_l{int(lamda * 100):05d}_m2m0400.bin"])
    end = time.time()
    print(f"Ellapsed time {end - start:.2f}s")
