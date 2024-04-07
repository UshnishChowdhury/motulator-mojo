from math import sqrt

alias PI = 3.141592653589793

@value
struct BaseValues:
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

    var U_nom: Float32
    var I_nom: Float32
    var f_nom: Float32
    var P_nom: Float32
    var tau_nom: Float32
    var n_p: Int
