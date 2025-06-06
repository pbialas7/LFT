#!/usr/bin/env python3

from concurrent.futures import ThreadPoolExecutor
import subprocess
import argparse
import math
import re

import numpy as np

beta_cr = math.log(1 + math.sqrt(2)) / 2


def beta_cnv(str):
    return beta_cr if str == 'cr' else float(str)


def eval_list(str, cnv=lambda t: float(t)):
    tok = str[1:-1].split(',')
    return [cnv(t) for t in tok]


def eval_range(str, cnv=lambda t: float(t)):
    m = range_reg.match(str)
    if m is None:
        raise ValueError(f"Not matching range format: {str}")
    tok = m.group(1).split(',')
    if len(tok) != 3:
        raise ValueError(f"Invalid range format: {str}")

    return np.arange(*[cnv(t) for t in tok])


list_reg = re.compile(r'\[(.*)\]')
range_reg = re.compile(r'range\((.*)\)')


def process_beta(str):
    if range_reg.match(str):
        return eval_range(str, cnv=beta_cnv)
    elif list_reg.match(str):
        return eval_list(str, cnv=beta_cnv)
    else:
        return [beta_cnv(str)]


def run_job(beta, args):
    subprocess.run([executable, "--beta", str(beta), '--Lx', str(args.L), '--Ly', str(args.L),
                    '--name', name, '-t', str(args.t), '-n', str(args.n), '--wolff',
                    '--data-dir', data_dir,
                    '--meas-freq', str(1), '--save-freq', str(0), '--seed', str(args.seed),
                    '--corr-freq',
                    str(args.corr_freq)])


executable = "../cmake-build-release/ising"
data_dir = "data"
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run ising model")
    parser.add_argument('-b', '--beta', type=str, help="Inverse temperature (beta)", required=True)
    parser.add_argument('-L', type=int, help="Lattice size (L)", required=True)
    parser.add_argument('-t', type=int, help="Number of Monte Carlo steps tuning sweeps", default=0)
    parser.add_argument('-n', type=int, help="Number of Monte Carlo sweeps measurement sweeps", default=0)
    parser.add_argument('--heath-bath', action='store_true', help="Use heath bath algorithm")
    parser.add_argument('--seed', type=int, help="Random seed", default=54746453)
    parser.add_argument('--corr-freq', type=int, help="Correlation frequency", default=10)
    args = parser.parse_args()

    beta_arg = process_beta(args.beta)

    print(beta_arg)
    with ThreadPoolExecutor(max_workers=4) as executor:
        for beta in beta_arg:
            beta_string = f"{round(beta * 1.e9):010d}"

            name = f"{args.L:04d}x{args.L:04d}_b{beta_string}"

            if args.heath_bath:
                name += "_hb"

            print(name)
            executor.submit(run_job, beta, args)
