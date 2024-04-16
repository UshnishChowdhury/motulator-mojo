
@value
struct ModelPars:
    """
    Inverse-Γ model parameters of an induction machine.
    
    Parameters
    ----------
    R_s : float
        Stator resistance (Ω).
    R_R : float
        Rotor resistance (Ω).
    L_sgm : float
        Leakage inductance (H).
    L_M : float
        Magnetizing inductance (H).
    n_p : int
        Number of pole pairs.  
    J : float
        Moment of inertia (kgm²).  
    
    """
    var R_s: Float16
    var R_R: Float16
    var L_sgm: Float16
    var L_M: Float16
    #var n_p: Int
    #var J: Float16