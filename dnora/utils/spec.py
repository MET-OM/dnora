import numpy as np
from scipy import interpolate


def interp_spec(f, D, S, fi, Di):
    """Interpolates a spectrum to new frequncy and directions.
    Spectrum is two dimensional (len(f), len(D))
    """
    Sleft = S
    Sright = S
    Dleft = -D[::-1]
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
    ind_flip = ((ind - 2 * steps * direction).astype(int) + len(D)) % len(D)

    spec_flip = spec[..., list(ind_flip)]

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

    ind_flip = ((ind + int(shift / dD)).astype(int) + len(D)) % len(D)

    spec_shift = spec[..., list(ind_flip)]
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
