from complex import ComplexSIMD
from python import Python


def main():
    var a  = 0
    var b  = 2

    var number1: ComplexSIMD[DType.float16, 1]
    number1 = number1.__init__(a, b)
    var number2: ComplexSIMD[DType.float16, 1]
    number2 = number2.__init__(a, b)

    var number3 = number1.__mul__(number2)

    print(number3.im)
    