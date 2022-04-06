
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
    res *= delta(j1, j2, j3)
    sumTerm = 0
    for i in range(min(j3-m3, j2-j1+j3) + 1):
        sumTerm += vterm(j1, j2, j3, m1, m2, m3, i)
    return res * sumTerm

def wigner6j(e, a, f, b, d, c):
    pre = delta(a,b,c)*delta(a,e,f)*delta(c,d,e)*delta(b,d,f)*fac(a+b+c+1)*fac(b+d+f+1)/(fac(a+b-c)*fac(c-d+e)*fac(c+d-e)*fac(c+d-e)*fac(a-e+f)*fac(-a+e+f)*fac(b+d-f))

    sumPart = 0
    for z in range(min([2*b,-a+b+c, b-d+f]) + 1):
        tmp = fac(2*b-z)*fac(b+c-e+f-z)*fac(b+c+e+f+1-z) / (fac(z)*fac(-a+b+c-z)*fac(b-d+f-z)*fac(a+b+c+1-z)*fac(b+d+f+1-z))
        sumPart += (1 if z%2==0 else -1) * tmp

    return (1 if (b+c+e+f)%2==0 else -1) * pre * sumPart

def cg(j1, j2, j3, m1, m2, m3):
    phase = 1 if int(-j1+j2-m3) % 2 == 0 else -1
    return phase*sqrt(2*j3+1)*wigner3j(j1, j2, j3, m1, m2, m3)

if __name__ == '__main__':
    print("Hello")

    # c = cg(1, 1, 2, 0, 0, 0)
    # print(c*c)

    # print(wigner6j(2, 1, 1, 2, 1, 1) * 30)
    print(wigner6j(9, 15, 21, 21, 21, 21))
