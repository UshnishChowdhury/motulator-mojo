from complex import ComplexSIMD
from python import Python

def main():
    var i_rs: ComplexSIMD[DType.float16, 1]
    i_rs.re = 1
    i_rs.im = 2
    print(i_rs.im)
    var np = Python.import_module("numpy")
    print(np.conj(i_rs.im))