from math import sqrt
from tensor import Tensor
from complex import ComplexSIMD
from python import Python

alias PI = 3.141592653589793


fn abc2complex(u: Tensor[DType.float16]) -> ComplexSIMD[DType.float16, 1]:
    """
    Transform three-phase quantities to a complex space vector.

    Parameters
    ----------
    u : array_like, shape (3,)
        Phase quantities.

    Returns
    -------
    complex
        Complex space vector (peak-value scaling).

    Examples
    --------
    >>> from motulator import abc2complex
    >>> y = abc2complex([1, 2, 3])
    >>> y
    (-1-0.5773502691896258j)

    """
    return convert_to_complex_form((2 / 3) * u[0] - (u[1] + u[2]) / 3).__add__(
        j().__mul__(convert_to_complex_form((u[1] - u[2]) / sqrt(3)))
    )


# %%
fn complex2abc(u: ComplexSIMD[DType.float16, 1]) raises -> PythonObject:
    """
    Transform a complex space vector to three-phase quantities.

    Parameters
    ----------
    u : complex
        Complex space vector (peak-value scaling).

    Returns
    -------
    ndarray, shape (3,)
        Phase quantities.

    Examples
    --------
    >>> from motulator import complex2abc
    >>> y = complex2abc(1-.5j)
    >>> y
    array([ 1.       , -0.9330127, -0.0669873])

    """

    var np = Python.import_module("numpy")
    return np.array(
        [u.re, 0.5 * (-u.re + sqrt(3) * u.im), 0.5 * (-u.re - sqrt(3) * u.im)]
    )


fn convert_to_complex_form(number: Float16) -> ComplexSIMD[DType.float16, 1]:
    var converted_complex_number: ComplexSIMD[DType.float16, 1]
    converted_complex_number = converted_complex_number.__init__(number, 0)
    return converted_complex_number


fn j() -> ComplexSIMD[DType.float16, 1]:
    var imag_j: ComplexSIMD[DType.float16, 1]
    imag_j = imag_j.__init__(0, 1)
    return imag_j


@value
struct NominalValues:
    """
    Base values.

    Base values are computed from the nominal values and the number of pole
    pairs. They can be used, e.g., for scaling the plotted waveforms.

    Parameters
    ----------
    U_nom : float
        Voltage (V, rms, line-line).
    I_nom : float
        Current (A, rms).
    f_nom : float
        Frequency (Hz).
    tau_nom : float
        Torque (Nm).
    P_nom : float
        Power (W).
    n_p : int
        Number of pole pairs.

    """

    var U_nom: Float32
    var I_nom: Float32
    var f_nom: Float32
    var P_nom: Float32
    var tau_nom: Float32
    var n_p: Int

    fn calculateBaseValuesFromNominalValues(self) -> BaseValues:
        var u = sqrt(2 / 3) * self.U_nom
        var i = sqrt(2) * self.I_nom
        var w = 2 * PI * self.f_nom
        var psi = u / w
        var p = 1.5 * u * i
        var Z = u / i
        var L = Z / w
        var tau = self.n_p * p / w

        return BaseValues(u, i, w, psi, p, Z, L, tau)


@value
struct BaseValues:
    """
    Attributes
    ----------
    u : float
        Base voltage (V, peak, line-neutral).
    i : float
        Base current (A, peak).
    w : float
        Base angular frequency (rad/s).
    psi : float
        Base flux linkage (Vs).
    p : float
        Base power (W).
    Z : float
        Base impedance (Ω).
    L : float
        Base inductance (H).
    tau : float
        Base torque (Nm).
    """

    var u: Float32
    var i: Float32
    var w: Float32
    var psi: Float32
    var p: Float32
    var Z: Float32
    var L: Float32
    var tau: Float32
