import torch

def plaquettes_batch(J: torch.Tensor) -> torch.Tensor:
    """
    Compute plaquettes for a batch of 2D lattices without using inverses.

    Args:
        J: Tensor of shape (B, 2, L, L)
           J[:, 0] -> x-links, J[:, 1] -> y-links

    Returns:
        Tensor of shape (B, L, L) with plaquette products:
        P(x,y) = Ux(x,y) * Uy(x+1,y) * Ux(x,y+1) * Uy(x,y)
        using periodic boundary conditions.
    """
    if J.ndim != 4 or J.shape[1] != 2 or J.shape[2] != J.shape[3]:
        raise ValueError(f"Expected J with shape (B, 2, L, L), got {tuple(J.shape)}")

    Ux = J[:, 0]  # (B, L, L)
    Uy = J[:, 1]  # (B, L, L)

    Uy_xplus = torch.roll(Uy, shifts=-1, dims=1)  # Uy(x+1, y)
    Ux_yplus = torch.roll(Ux, shifts=-1, dims=2)  # Ux(x, y+1)

    return Ux * Uy_xplus * Ux_yplus * Uy

import torch

def plaquettes_2x1_2d(J: torch.Tensor) -> torch.Tensor:
    """
    Compute 2x1 rectangular plaquettes on a 2D periodic lattice (no inverse).

    Args:
        J: Tensor of shape (B, 2, L, L)
           J[:, 0] = x-links, J[:, 1] = y-links

    Returns:
        Tensor of shape (B, L, L) with
        P_2x1(x,y) =
            Ux(x,y) * Ux(x+1,y) *
            Uy(x+2,y) *
            Ux(x+1,y+1) * Ux(x,y+1) *
            Uy(x,y)
    """
    if J.ndim != 4 or J.shape[1] != 2 or J.shape[2] != J.shape[3]:
        raise ValueError(f"Expected J with shape (B, 2, L, L), got {tuple(J.shape)}")

    Ux = J[:, 0]  # (B, L, L)
    Uy = J[:, 1]  # (B, L, L)

    # periodic shifts: lattice axes are dims 1 (x) and 2 (y)
    Ux_xp1 = torch.roll(Ux, shifts=-1, dims=1)      # Ux(x+1, y)
    Ux_yp1 = torch.roll(Ux, shifts=-1, dims=2)      # Ux(x, y+1)
    Ux_xp1_yp1 = torch.roll(Ux_yp1, shifts=-1, dims=1)  # Ux(x+1, y+1)
    Uy_xp2 = torch.roll(Uy, shifts=-2, dims=1)      # Uy(x+2, y)

    return Ux * Ux_xp1 * Uy_xp2 * Ux_xp1_yp1 * Ux_yp1 * Uy

def plaquettes_1x2_2d(J: torch.Tensor) -> torch.Tensor:
    """
    Compute 1x2 rectangular plaquettes on a 2D periodic lattice (no inverse).

    Args:
        J: Tensor of shape (B, 2, L, L)
           J[:, 0] = x-links, J[:, 1] = y-links

    Returns:
        Tensor of shape (B, L, L) with
        P_1x2(x,y) =
            Ux(x,y) *
            Uy(x+1,y) * Uy(x+1,y+1) *
            Ux(x,y+2) *
            Uy(x,y+1) * Uy(x,y)
    """
    if J.ndim != 4 or J.shape[1] != 2 or J.shape[2] != J.shape[3]:
        raise ValueError(f"Expected J with shape (B, 2, L, L), got {tuple(J.shape)}")

    Ux = J[:, 0]  # (B, L, L)
    Uy = J[:, 1]  # (B, L, L)

    # periodic shifts: lattice axes are dims 1 (x) and 2 (y)
    Uy_xp1 = torch.roll(Uy, shifts=-1, dims=1)          # Uy(x+1, y)
    Uy_yp1 = torch.roll(Uy, shifts=-1, dims=2)          # Uy(x, y+1)
    Uy_xp1_yp1 = torch.roll(Uy_xp1, shifts=-1, dims=2)  # Uy(x+1, y+1)
    Ux_yp2 = torch.roll(Ux, shifts=-2, dims=2)          # Ux(x, y+2)

    return Ux * Uy_xp1 * Uy_xp1_yp1 * Ux_yp2 * Uy_yp1 * Uy


def plaquettes_plane_3d(J: torch.Tensor, mu: int, nu: int) -> torch.Tensor:
    """
    Plaquettes in a single (mu, nu) plane for a batch of 3D lattices (no inverse).

    Args:
        J: Tensor of shape (B, 3, L, L, L)
           J[:, d] are links along direction d in {0,1,2}.
        mu, nu: Plane directions (different integers in {0,1,2}).

    Returns:
        Tensor of shape (B, L, L, L):
        P_{mu,nu}(x) = U_mu(x) * U_nu(x+mu) * U_mu(x+nu) * U_nu(x)
        with periodic boundary conditions.
    """
    if J.ndim != 5 or J.shape[1] != 3 or not (J.shape[2] == J.shape[3] == J.shape[4]):
        raise ValueError(f"Expected J with shape (B, 3, L, L, L), got {tuple(J.shape)}")
    if mu == nu or mu not in (0, 1, 2) or nu not in (0, 1, 2):
        raise ValueError(f"mu and nu must be different and in {{0,1,2}}, got mu={mu}, nu={nu}")

    # lattice dims in J are at axes 2,3,4; after selecting direction they become 1,2,3
    shift_dim_mu = 1 + mu
    shift_dim_nu = 1 + nu

    U_mu = J[:, mu]  # (B, L, L, L)
    U_nu = J[:, nu]  # (B, L, L, L)

    U_nu_xplusmu = torch.roll(U_nu, shifts=-1, dims=shift_dim_mu)  # U_nu(x+mu)
    U_mu_xplusnu = torch.roll(U_mu, shifts=-1, dims=shift_dim_nu)  # U_mu(x+nu)

    return U_mu * U_nu_xplusmu * U_mu_xplusnu * U_nu  