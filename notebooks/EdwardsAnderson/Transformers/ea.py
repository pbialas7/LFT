import torch


def Energy(spins, J):
    e = torch.sum(
        spins * torch.roll(spins, -1, -2) * J[0]
        , dim=(-2, -1)
    )

    e += torch.sum(
        spins * torch.roll(spins, -1, -1) * J[1]
        , dim=(-2, -1)
    )

    return -e
