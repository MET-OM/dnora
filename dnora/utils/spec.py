import numpy as np
from scipy import interpolate
import scipy


def directional_distribution(freq, fp, dirs, dirp):
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

    if np.max(dirs) < 7:
        raise Warning(
            f"Max value of directional vector is {np.max(dirs)}. Are you sure it is given in degrees?"
        )

    dirs = dirs - dirp

    mask = dirs > 180
    dirs[mask] = dirs[mask] - 360
    mask = dirs < -180
    dirs[mask] = dirs[mask] + 360

    theta = np.deg2rad(dirs)

    D = np.zeros((len(freq), len(theta)))

    for n in range(len(freq)):
        D[n, :] = A2[n] * np.cos(0.5 * theta) ** (2 * s[n])
    return D


def add_directional_distribution_to_spectrum(spec1d, D):
    """Add a directional distribution to the spectrum"""
    spec2d = np.zeros(D.shape)
    for n in range(D.shape[-1]):
        spec2d[:, n] = D[:, n] * spec1d

    return spec2d


def expand_to_directional_spectrum(spec1d, freq, dirs, dirp=None):
    inds = np.argmax(spec1d, axis=2)
    fp = freq[inds]
    if dirp is None:
        dirp = np.zeros(fp.shape)
    Nt, N, Nf = spec1d.shape
    spec2d = np.zeros(((Nt, N, Nf, len(dirs))))
    for n in range(N):
        for t in range(Nt):
            D = directional_distribution(freq, fp[t, n], dirs, dirp[t, n])
            spec2d[t, n, :, :] = add_directional_distribution_to_spectrum(
                spec1d[t, n, :], D
            )

    return spec2d


def _jonswap_one_spec(m0, fp, freq, gamma) -> np.ndarray:
    """Calculate one JONSWAP spectrum"""
    alpha = 1  # We will normalize anyway
    sigma_a = 0.07
    sigma_b = 0.09
    g = 9.81
    E_JS = (
        alpha * g**2 * (2 * np.pi) ** -4 * freq**-5 * np.exp(-5 / 4 * (freq / fp) ** -4)
    )
    sigma = np.full(len(freq), sigma_a)
    sigma[freq > fp] = sigma_b
    G_exp = np.exp(-0.5 * ((freq / fp - 1) / sigma) ** 2)
    E_JS = E_JS * gamma**G_exp
    var = np.trapz(E_JS, freq)
    E_JS = E_JS * m0 / var

    return E_JS


def jonswap1d(m0, fp, freq, gamma=3.3) -> np.ndarray:

    Nt, N = m0.shape
    E_JS = np.zeros((Nt, N, len(freq)))
    for n in range(N):
        for t in range(Nt):
            E_JS[t, n, :] = _jonswap_one_spec(m0[t, n], fp[t, n], freq, gamma)

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
    possible_shapes = []
    if spec.shape[-1] == len(dirs) and spec.shape[-2] == len(freq):
        possible_shapes.append(2)
    if spec.shape[-1] == len(freq) and spec.shape == dirs.shape:
        possible_shapes.append(1)

    if expected_dim is None:
        if 1 not in possible_shapes and 2 not in possible_shapes:
            raise ValueError("Provided array does not contain valid 1D or 2D spectra!")
        return possible_shapes

    if expected_dim not in possible_shapes:
        if 1 not in possible_shapes and 2 not in possible_shapes:
            raise ValueError("Provided array does not contain valid 1D or 2D spectra!")
        else:
            raise ValueError(
                f"Expected {expected_dim} dimensional spectra, but they seem to be {possible_shapes} dimensional!"
            )
