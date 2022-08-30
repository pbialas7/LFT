import numpy as np


def block_mean(sample, block_size):
    n_blocks = len(sample) // block_size
    if n_blocks == 0:
        raise RuntimeError("sample shorter then block_size")
    sample = sample.reshape(n_blocks, block_size, *sample.shape[1:])
    return sample.mean(1)


def ac(series, history=100):
    mu = np.mean(series)
    var = np.mean((series - mu) * (series - mu))
    out = np.empty(history)
    for i in range(history):
        out[i] = np.mean((series[:-history] - mu) * (series[i : i - history] - mu))
    return out / var


def ac_and_tau_int(series, c=10, maxlen=200):
    mu = np.mean(series)
    var = np.mean((series - mu) * (series - mu))
    out = [1.0]
    tau_int = 0.5
    for t in range(1, maxlen):
        cor = np.mean((series[:-t] - mu) * (series[t:] - mu)) / var
        tau_int += cor
        out.append(cor)
        if t > c * tau_int:
            break
    return tau_int, np.asarray(out)


def list_mean(data, mean=0.0):
    mean = mean
    size = 0
    for t in data:
        mean += t.sum()
        size += t.nelement()

    return mean / size
