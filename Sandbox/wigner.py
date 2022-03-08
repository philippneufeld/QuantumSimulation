
from numpy import sqrt

def fac(n):
    res = 1
    for i in range(n):
        res *= (i+1)
    return res

def delta(j1, j2, j3):
    return sqrt(fac(j1+j2-j3)*fac(j1-j2+j3)*fac(-j1+j2+j3)/fac(j1+j2+j3+1))

def vterm(j1, j2, j3, m1, m2, m3, v):
    ex = int(2*j2-j1-m1+v)
    phase = 1 if ex % 2 == 0 else -1
    u = j2-j1+j3
    return phase*(fac(j2+j3-m2-m3-v)*fac(j1+m2+m3+v)/(fac(v)*fac(j3-m3-v)*fac(u-v)*fac(j3+m3-u+v)))

def wigner3j(j1, j2, j3, m1, m2, m3):
    res = sqrt(fac(j3-m3)*fac(j3+m3)/(fac(j2+m2)*fac(j2-m2)*fac(j1-m2-m3)*fac(j1+m2+m3)))
    
    print(res)
    
    res *= delta(j1, j2, j3)
    
    print(delta(j1, j2, j3))
    print(res)

    sumTerm = 0
    for i in range(min(j3-m3, j2-j1+j3) + 1):
        sumTerm += vterm(j1, j2, j3, m1, m2, m3, i)

    print(sumTerm)

    return res * sumTerm

def cg(j1, j2, j3, m1, m2, m3):
    phase = 1 if int(-j1+j2-m3) % 2 == 0 else -1
    return phase*sqrt(2*j3+1)*wigner3j(j1, j2, j3, m1, m2, m3)

if __name__ == '__main__':
    print("Hello")

    c = cg(1, 1, 2, 0, 0, 0)
    print(c*c)