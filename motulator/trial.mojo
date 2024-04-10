from complex import ComplexSIMD

def main():
    var i_rs: ComplexSIMD[DType.float16, 1]
    i_rs.re = 1
    i_rs.im = 2
    print(i_rs.im)