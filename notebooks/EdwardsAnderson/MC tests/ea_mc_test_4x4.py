# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.18.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
# %load_ext autoreload
# %autoreload 2

# %%
import numpy as np
import torch
import matplotlib.pyplot as plt

# %%
import yaml

# %%
import dneumc

# %%
from dneumc.ising.ising import fast_gen_states

# %%
from dneumc.ising.ea import Energy

# %%
from neumc.utils.stats_utils import bootstrap, ac_and_tau_int, torch_bootstrap

# %%
ea_dir = '../../../ea_pt_mt_test'

# %%
Lx = 4
Ly = 6

# %%
spins = fast_gen_states(Lx*Ly).view(-1,Lx,Ly)

# %%
n_replicas = 2

# %%
name =f'{Lx:02d}x{Ly:02d}_bxxxx_ex01'

# %%
with open(f"{ea_dir}/opt_{name}.yaml", "r") as f:
    options = yaml.safe_load(f)
n_chains=n_betas=len(options['beta'])
betas = options['beta']
print(betas)

# %%
J = torch.from_numpy(np.loadtxt(f"{ea_dir}/j_{name}.txt").reshape(2,Lx,Ly))
J

# %%
energy = Energy(spins,J)

# %%
for i,beta in enumerate(betas):
    print(f'beta = {beta}') 

    cfgs = torch.from_numpy(np.fromfile(f"{ea_dir}/cfg_{name}_b{i:02d}.bin", dtype="int8").reshape(-1,n_replicas,Lx,Ly))
    em = torch.from_numpy(np.loadtxt(f"{ea_dir}/em_{name}_b{i:02d}.txt")) 

    n_cfg = min(100000,cfgs.shape[0])
    for r in range(n_replicas):
        if not torch.allclose(em[:n_cfg,2*r],Energy(cfgs[:n_cfg,r],J)):
            print(f'Energies for beta = {beta} in replica {r} do not match')

    print(f'tau ',end='')
    for r in range(n_replicas):
        tau, ac = ac_and_tau_int(em[:,2*r+1].numpy());
        print(f'M[{r}] {tau:.2f} ', end='')
    print()    
    
    Z = torch.exp(-beta*energy).sum()
    e_th = torch.sum(energy*torch.exp(-beta*energy))/Z
    print('E ',end='')
    for r in range(n_replicas):
        e_val, e_err = torch_bootstrap(em[:,2*r],n_samples=200, binsize=10)
        print(f"{e_val-e_th:.4f} {e_err:.4f} ",end='')
        if np.abs(e_val.item()-e_th.item())/e_err.item() >3 :
            print('3 sigma disrepancy in energy')
    print()    

    mag = spins.sum((-1,-2))
    am_th = torch.sum(torch.abs(mag)*torch.exp(-beta*energy))/Z
    print('M ',end='')
    for r in range(n_replicas):
        am_val, am_err = torch_bootstrap(torch.abs(em[:,2*r+1]),n_samples=200, binsize=10)
        print(f"{am_val-am_th:.4f} {am_err:.4f} ",end='')
        if np.abs(am_val.item()-am_th.item())/am_err.item() >3 :
            print('3 sigma disrepancy in |M|')
    print()  
        
        

# %%

# %%
