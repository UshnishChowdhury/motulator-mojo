from complex import ComplexSIMD


struct InductionMachine:
    """
    Γ-equivalent model of an induction machine.

    An induction machine is modeled using the Γ-equivalent model [#Sle1989]_.
    The model is implemented in stator coordinates. The flux linkages are used
    as state variables.

    Parameters
    ----------
    n_p : int
        Number of pole pairs.
    R_s : float
        Stator resistance (Ω).
    R_r : float
        Rotor resistance (Ω).
    L_ell : float
        Leakage inductance (H).
    L_s : float
        Stator inductance (H).

    Notes
    -----
    The Γ model is chosen here since it can be extended with the magnetic
    saturation model in a staightforward manner. If the magnetic saturation is
    omitted, the Γ model is mathematically identical to the inverse-Γ and T
    models [#Sle1989]_.

    References
    ----------
    .. [#Sle1989] Slemon, "Modelling of induction machines for electric drives,"
       IEEE Trans. Ind. Appl., 1989, https://doi.org/10.1109/28.44251.

    """

    var n_p: Int
    var R_s: Float16
    var R_r: Float16
    var L_ell: Float16
    var L_s: Float16
    var psi_ss0: ComplexSIMD[DType.float16, 1]
    var psi_rs0: ComplexSIMD[DType.float16, 1]

    fn __init__(
        inout self, n_p: Int, R_s: Float16, R_r: Float16, L_ell: Float16, L_s: Float16
    ):
        self.n_p = n_p
        self.R_s, self.R_r = R_s, R_r
        self.L_ell, self.L_s = L_ell, L_s
        # Initial values
        self.psi_ss0.im = 0
        self.psi_rs0.im = 0

    fn currents(
        inout self,
        psi_ss: ComplexSIMD[DType.float16, 1],
        psi_rs: ComplexSIMD[DType.float16, 1],
    ) -> ComplexCurrents:
        """
        Compute the stator and rotor currents.

        Parameters
        ----------
        psi_ss : complex
            Stator flux linkage (Vs).
        psi_rs : complex
            Rotor flux linkage (Vs).

        Returns
        -------
        i_ss : complex
            Stator current (A).
        i_rs : complex
            Rotor current (A).

        """
        var i_rs: ComplexSIMD[DType.float16, 1]
        var i_ss: ComplexSIMD[DType.float16, 1]

        i_rs.re = (psi_rs.re - psi_ss.re) / (self.L_ell)
        i_rs.im = (psi_rs.im - psi_ss.im) / (self.L_ell)

        i_ss.re = psi_ss.re / self.L_s - i_rs.re
        i_ss.im = psi_ss.im / self.L_s - i_rs.im

        return ComplexCurrents(i_ss, i_rs)

    fn magnetic(
        inout self,
        psi_ss: ComplexSIMD[DType.float16, 1],
        psi_rs: ComplexSIMD[DType.float16, 1],
    ) -> MagneticModel:
        """
        Magnetic model.

        Parameters
        ----------
        psi_ss : complex
            Stator flux linkage (Vs).
        psi_rs : complex
            Rotor flux linkage (Vs).

        Returns
        -------
        i_ss : complex
            Stator current (A).
        i_rs : complex
            Rotor current (A).
        tau_M : float
            Electromagnetic torque (Nm).

        """
        var complex_currents: ComplexCurrents = self.currents(psi_ss, psi_rs)
        var current_psi_ss = psi_ss
        # Changing imaginary axis sign for complex conjugate
        current_psi_ss.im = -1 * psi_ss.im
        var tau_M = 1.5 * self.n_p * (complex_currents.i_ss.__mul__(current_psi_ss)).im

        return MagneticModel(complex_currents.i_ss, complex_currents.i_rs, tau_M)

    fn f(
        inout self,
        psi_ss: ComplexSIMD[DType.float16, 1],
        psi_rs: ComplexSIMD[DType.float16, 1],
        u_ss: ComplexSIMD[DType.float16, 1],
        w_M: Float16,
    ):
        """
        Compute the state derivatives.

        Parameters
        ----------
        psi_ss : complex
            Stator flux linkage (Vs).
        psi_rs : complex
            Rotor flux linkage (Vs).
        u_ss : complex
            Stator voltage (V).
        w_M : float
            Rotor angular speed (mechanical rad/s).

        Returns
        -------
        complex list, length 2
            Time derivative of the state vector, [dpsi_ss, dpsi_rs]
        i_ss : complex
            Stator current (A).
        tau_M : float
            Electromagnetic torque (Nm).

        Notes
        -----
        In addition to the state derivatives, this method also returns the
        output signals (stator current `i_ss` and torque `tau_M`) needed for
        interconnection with other subsystems. This avoids overlapping
        computation in simulation.

        """
        # var i_ss, i_rs, tau_M = self.magnetic(psi_ss, psi_rs)
        var magnetic_model = self.magnetic(psi_ss, psi_rs)
        var dpsi_ss: ComplexSIMD[DType.float16, 1]
        dpsi_ss = dpsi_ss.__init__(
            self.calculate_dpsi(magnetic_model.i_ss.re, u_ss.re),
            self.calculate_dpsi(magnetic_model.i_ss.im, u_ss.im),
        )

        var dpsi_rs = (
            -self.convert_to_complex_form(self.R_r).__mul__(magnetic_model.i_rs)
        ).__add__(
            self.j().__mul__(
                self.convert_to_complex_form(self.n_p)
                .__mul__(self.convert_to_complex_form(w_M))
                .__mul__(psi_rs)
            )
        )

        return [dpsi_ss, dpsi_rs], magnetic_model.i_ss, magnetic_model.tau
    
    def meas_currents(self):
        """
        Measure the phase currents at the end of the sampling period.

        Returns
        -------
        i_s_abc : 3-tuple of floats
            Phase currents (A).

        """
        # Stator current space vector in stator coordinates
        var i_ss = self.currents(self.psi_ss0, self.psi_rs0)
        # Phase currents
        i_s_abc = complex2abc(i_ss)  # + noise + offset ...
        return i_s_abc

    fn calculate_dpsi(inout self, i_ss: Float16, u_ss: Float16) -> Float16:
        return u_ss - self.R_s * i_ss

    fn convert_to_complex_form(
        inout self, number: Float16
    ) -> ComplexSIMD[DType.float16, 1]:
        var converted_complex_number: ComplexSIMD[DType.float16, 1]
        converted_complex_number = converted_complex_number.__init__(number, 0)
        return converted_complex_number

    fn j(inout self) -> ComplexSIMD[DType.float16, 1]:
        var imag_j: ComplexSIMD[DType.float16, 1]
        imag_j = imag_j.__init__(0, 1)
        return imag_j


struct InductionMachineInvGamma[InductionMachine]:
    """
    Inverse-Γ model of an induction machine.

    This extends the InductionMachine class (based on the Γ model) by providing
    an interface for the inverse-Γ model parameters. Linear magnetics are
    assumed. If magnetic saturation is to be modeled, the Γ model is preferred.

    Parameters
    ----------
    n_p : int
        Number of pole pairs.
    R_s : float
        Stator resistance (Ω).
    R_R : float
        Rotor resistance (Ω).
    L_sgm : float
        Leakage inductance (H).
    L_M : float
        Magnetizing inductance (H).

    """

    def __init__(self, n_p, R_s, R_R, L_sgm, L_M):
        # pylint: disable=too-many-arguments
        # Convert the inverse-Γ parameters to the Γ parameters
        gamma = L_M / (L_M + L_sgm)  # Magnetic coupling factor
        super().__init__(n_p, R_s, R_R / gamma**2, L_sgm / gamma, L_M + L_sgm)
        # Initial values
        self.psi_ss0 = 0j
        self.psi_rs0 = 0j  # self.psi_rs0 = self.psi_Rs0/gamma


@value
struct ComplexCurrents:
    var i_rs: ComplexSIMD[DType.float16, 1]
    var i_ss: ComplexSIMD[DType.float16, 1]


@value
struct MagneticModel:
    var i_ss: ComplexSIMD[DType.float16, 1]
    var i_rs: ComplexSIMD[DType.float16, 1]
    var tau: Float16
