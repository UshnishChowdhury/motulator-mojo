struct Mechanics:
    """
    Mechanics subsystem.

    This models an equation of motion for stiff mechanics.

    Parameters
    ----------
    J : float
        Total moment of inertia (kgm²).
    tau_L_w : callable
        Load torque (Nm) as a function of speed, `tau_L_w(w_M)`. For example,
        ``tau_L_w = b*w_M``, where `b` is the viscous friction coefficient. The
        default is zero, ``lambda w_M: 0*w_M``.
    tau_L_t : callable
        Load torque (Nm) as a function of time, `tau_L_t(t)`. The default is 
        zero, ``lambda t: 0*t``.

    """

    def __init__(self, J, tau_L_w=lambda w_M: 0*w_M, tau_L_t=lambda t: 0*t):
        self.J = J
        self.tau_L_t, self.tau_L_w = tau_L_t, tau_L_w
        # Initial values
        self.w_M0, self.theta_M0 = 0, 0

    def f(self, t, w_M, tau_M):
        """
        Compute the state derivatives.

        Parameters
        ----------
        t : float
            Time (s).
        w_M : float
            Rotor angular speed (mechanical rad/s).
        tau_M : float
            Electromagnetic torque (Nm).

        Returns
        -------
        list, length 2
            Time derivatives of the state vector.

        """
        # Total load torque
        tau_L = self.tau_L_w(w_M) + self.tau_L_t(t)
        # Time derivatives
        dw_M = (tau_M - tau_L)/self.J
        dtheta_M = w_M
        return [dw_M, dtheta_M]

    def meas_speed(self):
        """
        Measure the rotor speed.

        This returns the rotor speed at the end of the sampling period.

        Returns
        -------
        w_M0 : float
            Rotor angular speed (mechanical rad/s).

        """
        # The quantization noise of an incremental encoder could be modeled
        # by means of the rotor angle, to be done later.
        return self.w_M0

    def meas_position(self):
        """
        Measure the rotor angle.

        This returns the rotor angle at the end of the sampling period.

        Returns
        -------
        theta_M0 : float
            Rotor angle (mechanical rad).

        """
        return self.theta_M0