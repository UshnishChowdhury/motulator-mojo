# %%


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

    var J: Float16
    var tau_L_t: Float16
    var tau_L_w: Float16
    var w_M0: Float16
    var theta_M0: Float16

    fn __init__(inout self, J: Float16):
        self.J = J
        # Initial values
        self.tau_L_t = 0
        self.tau_L_w = 0
        self.w_M0 = 0
        self.theta_M0 = 0

    fn f(
        inout self, t: Float16, w_M: Float16, tau_M: Float16
    ) -> ListLiteral[Float16, Float16]:
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
        var tau_L: Float16 = self.calculate_load_torque_as_fn_of_speed(
            0, w_M
        ) + self.calculate_load_torque_as_fn_of_time(0, t)
        # Time derivatives
        var dw_M: Float16 = (tau_M - tau_L) / self.J
        var dtheta_M: Float16 = w_M
        return [dw_M, dtheta_M]

    fn meas_speed(inout self) -> Float16:
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

    fn meas_position(self) -> Float16:
        """
        Measure the rotor angle.

        This returns the rotor angle at the end of the sampling period.

        Returns
        -------
        theta_M0 : float
            Rotor angle (mechanical rad).

        """
        return self.theta_M0

    fn calculate_load_torque_as_fn_of_speed(
        inout self, w_M: Float16, b: Float16
    ) -> Float16:
        return w_M * b

    fn calculate_load_torque_as_fn_of_time(
        inout self, t: Float16, b: Float16
    ) -> Float16:
        return t * b
