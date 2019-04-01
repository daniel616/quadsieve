from sieve import *


def gcd_test():
    assert gcd(73,3)==1
    assert gcd(20,10)==10
    assert gcd(5,5)==5
    assert gcd(-5,5)==5
    assert gcd( 144,56)==8

if __name__ == "__main__":
    gcd_test()
