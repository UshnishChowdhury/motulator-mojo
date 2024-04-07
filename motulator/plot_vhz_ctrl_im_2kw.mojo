from _helpers import BaseValues, Attributes

def main():
    var base = BaseValues(U_nom=400, I_nom=5, f_nom=50, tau_nom=14.6, P_nom=2.2e3, n_p=2)
    var attributes = BaseValues.calculateAttributeValues(base)
    print(attributes.u)