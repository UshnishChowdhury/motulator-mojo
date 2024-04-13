from _helpers import BaseValues, NominalValues
from _drive import InductionMachineInvGamma
from _mechanics import Mechanics

fn main():
    # %%
    # Setting the nominal values.
    var nominalValues = NominalValues(U_nom=400, I_nom=5, f_nom=50, tau_nom=14.6, P_nom=2.2e3, n_p=2)
    
    # Compute base values based on the nominal values.
    var baseValues = NominalValues.calculateBaseValuesFromNominalValues(nominalValues)

    # Machine model, using its inverse-Γ parameters
    var machine = InductionMachineInvGamma(R_s=3.7, R_R=2.1, L_sgm=.021, L_M=.224, n_p=2)

    # Mechanics model
    var mechanics = Mechanics(J=.015)
    
    print(mechanics.tau_L_t)