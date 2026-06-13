import numpy as np
import torch
import yaml
from scipy.stats import norm
from neumc.utils.stats_utils import ac_and_tau_int, torch_bootstrap


def check_norm(a, b):
    va, ea = a
    vb, eb = b
    std = np.sqrt(ea * ea + eb * eb)
    stat = (va - vb) / std
    p = 2 * norm.sf(np.abs(stat))
    return p, stat


def check_norm_v(x):
    a = (x[0], x[1])
    b = (x[2], x[3])
    return check_norm(a, b)


def read_options(data_dir, name):
    with open(f"{data_dir}/opt_{name}.yaml", "r") as f:
        options = yaml.safe_load(f)
    return options


def read_em_files(data_dir, name, betas):
    em = []
    for i, _ in enumerate(betas):
        em.append(torch.from_numpy(np.loadtxt(f"{data_dir}/em_{name}_b{i:02d}.txt")))
    lengths = [len(x) for x in em]
    min_length = min(lengths)
    em_ = [x[:min_length] for x in em]
    return torch.stack(em_, dim=0)


def get_stats(em, betas):
    n_betas = len(betas)
    n_obs = 3
    stats = np.empty((n_betas, 2, n_obs, 2))  # chains, replicas, observables, (mean, error)
    for i, beta in enumerate(betas):
        for k in (0, 1):
            # integrated autocorrelation time
            tau, _ = ac_and_tau_int(em[i][::, 2 * k + 1].numpy())
            stats[i, k, 0, 0] = tau
            stats[i, k, 0, 1] = 0

            # energy
            v, e = torch_bootstrap(em[i][::, 2 * k], n_samples=100, binsize=10)
            stats[i, k, 1, 0] = v.item()
            stats[i, k, 1, 1] = e.item()

            # absolute magnetization
            v, e = torch_bootstrap(torch.abs(em[i][::, 2 * k + 1]), n_samples=100, binsize=10)
            stats[i, k, 2, 0] = v.item()
            stats[i, k, 2, 1] = e.item()

    return stats


def compare(s1, s2):
    s = np.concat((s1, s2), axis=1)
    return np.apply_along_axis(check_norm_v, 1, s)


def check(p, *, threshold=0.01):
    suspicious = p < threshold
    if np.any(suspicious):
        print("Suspiciously different results:")
        print(p[suspicious])
        return False
    else:
        print("Results are consistent.")
        return True
