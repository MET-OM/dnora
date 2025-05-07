import numpy as np
from scipy import interpolate
import scipy


def directional_distribution(freq, fp, dirs):
    """Calculates directional cos**2s(0.5*theta) distribution of spectrum

    dirs given in degrees"""
    # Eq. 6.3.22 in Holthuisen (2007), base on observations from Young et al. (1996) and Ewans (1998)
    sigma = 26.9 * (freq / fp) ** 0.68
    mask = freq < fp
    sigma[mask] = 26.9 * (freq[mask] / fp) ** -1.05

    # Inversion of Eq. 6.3.26 in Holthuisen (2007)
    # Width parameter
    s = 2 / np.deg2rad(sigma) ** 2 - 1
    s = np.maximum(s, 1)
    A2 = scipy.special.gamma(s + 1) / (
        scipy.special.gamma(s + 1 / 2) * 2 * np.sqrt(np.pi)
    )
    theta = np.deg2rad(dirs)

    D = np.zeros((len(freq), len(theta)))

    for n in range(len(freq)):
        D[n, :] = A2[n] * np.cos(0.5 * theta) ** (2 * s[n])

    return D


def add_directional_distribution_to_spectrum(spec1d, D):
    """Add a directional distribution to the spectrum"""
    spec2d = np.zeros(D.shape)
    for n in range(len(theta)):
        spec2d[:, n] = D[:, n] * spec1d

    return spec2d


def jonswap1d(m0, fp, freq, gamma=3.3) -> np.ndarray:
    alpha = 1  # We will normalize anyway
    sigma_a = 0.07
    sigma_b = 0.09
    g = 9.81

    E_JS = np.zeros((len(m0), len(freq)))
    for n, (f, m) in enumerate(zip(fp, m0)):
        E_JS[n, :] = (
            alpha
            * g**2
            * (2 * np.pi) ** -4
            * freq**-5
            * np.exp(-5 / 4 * (freq / f) ** -4)
        )
        sigma = np.full(len(freq), sigma_a)
        sigma[freq > f] = sigma_b
        G_exp = np.exp(-0.5 * ((freq / f - 1) / sigma) ** 2)
        E_JS[n, :] = E_JS[n, :] * gamma**G_exp
        var = np.trapz(E_JS[n, :], freq)
        E_JS[n, :] = E_JS[n, :] * m / var

    return E_JS


def interp_spec(f, D, S, fi, Di):
    """Interpolates a spectrum to new frequncy and directions.
    Spectrum is two dimensional (len(f), len(D))
    """
    Sleft = S
    Sright = S
    Dleft = D - 360  # -D[::-1]
    Dright = D + 360

    bigS = np.concatenate((Sleft, S, Sright), axis=1)
    bigD = np.concatenate((Dleft, D, Dright))
    Finterpolator = interpolate.RectBivariateSpline(f, bigD, bigS, kx=1, ky=1, s=0)
    Si = Finterpolator(fi, Di)

    return Si


# -----------------------------------------------------------------------------


def flip_spec(spec, D):
    """Flips the directionality of the spectrum (clock/anticlockwise).
    To flip the directional vector, use flip_spec(D,D)
    """

    # This check enables us to flip directions with flip_spec(D,D)
    if len(spec.shape) == 1:
        flipping_dir = True
        spec = np.array([spec])
    else:
        flipping_dir = False
    spec_flip = np.zeros(spec.shape)

    ind = np.arange(0, len(D), dtype="int")
    # dD = np.diff(D).mean()
    dD = 360 / len(D)
    steps = D / dD  # How many delta-D from 0

    # Need to move indeces the other way if the vector is decreasing than if it is increasing
    direction = np.sign(np.median(np.diff(D)))

    ind_flip = (np.round((ind - 2 * steps * direction)) + len(D)) % len(D)

    spec_flip = spec[..., list(ind_flip.astype(int))]

    if flipping_dir:
        spec_flip = spec_flip[0]

    return spec_flip


def shift_spec(spec, D, shift=0):
    """Shifts the spectrum D degree. To shift the directional vector, use
    shift_spec(D, D, shift)
    """

    # This check enables us to flip directions with shift_spec(D, D, shift)
    if len(spec.shape) == 1:
        shifting_dir = True
        spec = np.array([spec])
    else:
        shifting_dir = False
    spec_shift = np.zeros(spec.shape)

    ind = np.arange(0, len(D), dtype="int")
    dD = 360 / len(D)

    if abs(np.floor(shift / dD) - shift / dD) > 0:
        raise Exception(
            f"Shift {shift} needs to be multiple of frequency resolution {dD}, but shift/dD={shift/dD}! Otherwise interpolation would be needed."
        )

    ind_flip = (np.round(ind + int(shift / dD)) + len(D)) % len(D)

    spec_shift = spec[..., list(ind_flip.astype(int))]
    if shifting_dir:
        spec_shift = spec_shift[0]

    return spec_shift


def check_that_spectra_are_consistent(
    spec, dirs, freq, expected_dim: int = None
) -> int:
    if spec.shape[-1] == len(dirs) and spec.shape[-2] == len(freq):
        spec_dim = 2
    elif spec.shape[-1] == len(freq) and spec.shape[-1] == len(dirs):
        spec_dim = 1
    else:
        spec_dim = -1

    if expected_dim is None:
        return spec_dim

    if spec_dim != expected_dim:
        if spec_dim not in [1, 2]:
            ValueError("Provided array does not contain valid 1D or 2D spectra!")
        else:
            ValueError(
                f"Expected {expected_dim} dimensional spectra, but they seem to be {spec_dim} dimensional!"
            )
