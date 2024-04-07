from _helpers import BaseValues, NominalValues

def main():
    # %%
    # Setting the nominal values.
    var nominalValues = NominalValues(U_nom=400, I_nom=5, f_nom=50, tau_nom=14.6, P_nom=2.2e3, n_p=2)
    
    # Compute base values based on the nominal values.
    var baseValues = NominalValues.calculateBaseValuesFromNominalValues(nominalValues)
    
    print(baseValues.u)