# %%
import argparse

from rich import print
from rich.console import Console

console = Console()

import numpy as np
from scipy.stats import norm

import torch


# %%
import yaml

# %%
import dneumc

# %%
from dneumc.ising.ising import fast_gen_states

# %%
from dneumc.ising.ea import Energy

# %%
from neumc.utils.stats_utils import ac_and_tau_int, torch_bootstrap

def compare(val, err, val_th):
    diff = val-val_th
    s = np.abs(diff)/err
    p = 2*norm.sf(s)
    if p <0.01:
        print(f'[red]{diff:.4f} {err:.4f} {p:.2e}[/red]|', end='')
    else:
        print(f'{diff:.4f} {err:.4f} {p:.2e}|', end='')



# %%
ea_dir = '../../../ea_pt_mt_test'

# %%

parser = argparse.ArgumentParser()
parser.add_argument('--Lx', type=int, required=True)
parser.add_argument('--Ly', type=int, required=True)
parser.add_argument('--n_replicas', type=int, default=2)
parser.add_argument('--tag', type=str, required=True)
args = parser.parse_args()

Lx = args.Lx
Ly = args.Ly
tag = args.tag

# %%
spins = fast_gen_states(Lx*Ly).view(-1,Lx,Ly)

# %%
n_replicas = 2

# %%
name =f'{Lx:02d}x{Ly:02d}_{tag}'

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
        compare(e_val.item(), e_err.item(), e_th.item())
    print()    

    mag = spins.sum((-1,-2))
    am_th = torch.sum(torch.abs(mag)*torch.exp(-beta*energy))/Z
    print('M ',end='')
    for r in range(n_replicas):
        am_val, am_err = torch_bootstrap(torch.abs(em[:,2*r+1]),n_samples=200, binsize=10)
        compare(am_val.item(), am_err.item(), am_th.item())
    print()  
        
        

# %%

# %%
