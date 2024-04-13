from complex import ComplexSIMD
from _helpers import convert_to_complex_form, get_complex_conjugate


@value
struct Inverter:
    """
    Inverter with constant DC-bus voltage and switching-cycle averaging.

    Parameters
    ----------
    u_dc : float
        DC-bus voltage (V).

    """

    var u_dc: Float16

    @staticmethod
    fn ac_voltage(
        q: ComplexSIMD[DType.float16, 1], u_dc: Float16
    ) -> ComplexSIMD[DType.float16, 1]:
        """
        Compute the AC-side voltage of a lossless inverter.

        Parameters
        ----------
        q : complex
            Switching state vector.
        u_dc : float
            DC-bus voltage (V).

        Returns
        -------
        u_ac : complex
            AC-side voltage (V).

        """
        var u_ac = q.__mul__(convert_to_complex_form(u_dc))
        return u_ac

    @staticmethod
    fn dc_current(
        q: ComplexSIMD[DType.float16, 1], i_ac: ComplexSIMD[DType.float16, 1]
    ) -> Float16:
        """
        Compute the DC-side current of a lossless inverter.

        Parameters
        ----------
        q : complex
            Switching state vector.
        i_ac : complex
            AC-side current (A).

        Returns
        -------
        i_dc : float
            DC-side current (A).

        """
        var i_dc = 1.5 * (q.__mul__(get_complex_conjugate(i_ac))).re
        return i_dc

    fn meas_dc_voltage(inout self) -> Float16:
        """
        Measure the DC-bus voltage.

        Returns
        -------
        float
            DC-bus voltage (V).

        """
        return self.u_dc
