# Philipp Neufeld, 2021-2022

import numpy as np
import sympy as sp
import scipy.integrate
import matplotlib.pyplot as plt

x0 = 0

def func(x):
    return 1 / ((x-x0)**2 + 1) / (np.pi)


def adaptive1Helper(xs, f, a, dx, depth, f0, f2, f4, atol, rtol, maxdepth):
    I1 = (dx/6) * (f0 + 4*f2 + f4)
    fevs = 2

    f1, f3 = f(a + 0.25*dx), f(a + 0.75*dx)
    I2 = (dx/12)*(f0+4*(f1+f3)+2*f2+f4)
    Ierr = (I2-I1) / 15
    I = I2 + Ierr # equivalent to Boole's rule
    
    if (depth > 1):
        xs.append(a + 0.25*dx)
        xs.append(a + 0.75*dx)

    if (np.abs(Ierr)-atol > rtol*np.abs(I) and depth < maxdepth):
        dx /= 2
        Isub1, fevs1 = adaptive1Helper(xs, f, a, dx, depth+1, f0, f1, f2, atol/2, rtol, maxdepth)
        Isub2, fevs2 = adaptive1Helper(xs, f, a+dx, dx, depth+1, f2, f3, f4, atol/2, rtol, maxdepth)
        I = Isub1 + Isub2
        fevs += fevs1 + fevs2

    return I, fevs

def adaptive1(f, a, b, n, atol, rtol, maxdepth):

    n -= (n-1) % 4
    secs = (n-1) // 4
    dx = (b-a)/secs

    atol /= secs

    f0, f1, f2 = f(a), f(a+0.5*dx), f(a+dx)
    xs = []

    I, fevs = adaptive1Helper(xs, f, a, dx, 1, f0, f1, f2, atol, rtol, maxdepth)
    fevs += 3

    for i in range(1, secs):
        f0, f1, f2 = f2, f(a+(i+0.5)*dx), f(a+(i+1)*dx)
        subI, subFevs = adaptive1Helper(xs, f, a+i*dx, dx, 1, f0, f1, f2, atol, rtol, maxdepth)
        I += subI
        fevs += 2 + subFevs

    return I, fevs, np.array(xs)


def adaptive2Helper(Iest, xs, f, a, dx, depth, f0, f2, f4, atol, rtol, maxdepth):
    I1 = (dx/6) * (f0 + 4*f2 + f4)
    fevs = 2

    f1, f3 = f(a + 0.25*dx), f(a + 0.75*dx)
    I2 = (dx/12)*(f0+4*(f1+f3)+2*f2+f4)
    Ierr = (I2-I1) / 15
    I = I2 + Ierr # equivalent to Boole's rule

    if (depth > 0):
        xs.append(a + 0.25*dx)
        xs.append(a + 0.75*dx)

    # Iest = np.abs(I)
    Iest = np.max([np.abs(I), np.abs(Iest)])
    if (np.abs(Ierr) - atol > rtol*Iest and depth < maxdepth):
    # if (np.abs(Ierr)-atol > rtol*np.abs(I)):
        dx /= 2
        Isub1, fevs1 = adaptive2Helper(Iest / 2, xs, f, a, dx, depth+1, f0, f1, f2, atol/2, rtol, maxdepth)
        Isub2, fevs2 = adaptive2Helper(Iest / 2, xs, f, a+dx, dx, depth+1, f2, f3, f4, atol/2, rtol, maxdepth)
        I = Isub1 + Isub2
        fevs += fevs1 + fevs2

    return I, fevs

def adaptive2(f, a, b, n, atol, rtol, maxdepth):

    n -= (n-1) % 4
    secs = (n-1) // 4
    dx = (b-a)/secs

    atol /= secs

    fs = f(np.array([a+0.5*i*dx for i in range(1+2*secs)]))
    Iest = abs((0.5*(fs[0]+fs[-1]) + np.sum(fs[1:-1]))*0.5*dx / secs)
    fevs = 1+2*secs

    xs = [a+0.5*i*dx for i in range(1+2*secs)]
    I, subFevs = adaptive2Helper(Iest, xs, f, a, dx, 1, fs[0], fs[1], fs[2], atol, rtol, maxdepth)
    fevs += subFevs

    for i in range(1, secs):
        subI, subFevs = adaptive2Helper(Iest, xs, f, a+i*dx, dx, 1, fs[2*i], fs[2*i+1], fs[2*i+2], atol, rtol, maxdepth)
        I += subI
        fevs += subFevs

    return I, fevs, np.array(xs)


if __name__ == "__main__":

    a = -250
    b = 250
    atol = 1e-5
    rtol = 1e-4
    depth = 100
    n = 40

    # calculate exact solution
    xsp = sp.Symbol("x")
    exact = float(sp.integrate(1/((xsp-x0)**2+1)/(sp.pi), (xsp, a, b)))
    print(f"Exact: {exact}")

    # adaptive 1
    # print("")
    # I1, fevs1, xs1 = adaptive1(func, a, b, n, atol, rtol, depth)
    # print(f"Adaptive1 {abs(I1 - exact)/exact:e} ({fevs1})")
    # x1 = np.linspace(a, b, fevs1)
    # print(f"Euler ({fevs1}): {abs(np.sum(func(x1)) * (b-a)/(fevs1-1) - exact)/exact:e}")
    # print(f"Simpson ({fevs1}): {abs(scipy.integrate.simpson(func(x1), x1) - exact)/exact:e}")

    # adaptive 2
    print("")
    I2, fevs2, xs2 = adaptive2(func, a, b, n, atol, rtol, depth)
    print(f"Adaptive2 {abs(I2 - exact)/exact:e} ({fevs2})")
    x2 = np.linspace(a, b, fevs2)
    print(f"Euler ({fevs2}): {abs(np.sum(func(x2)) * (b-a)/(fevs2-1) - exact)/exact:e}")
    print(f"Simpson ({fevs2}): {abs(scipy.integrate.simpson(func(x2), x2) - exact)/exact:e}")

    x = np.linspace(a, b, 10000)

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, gridspec_kw={'height_ratios': [4 ,1, 1]}, sharex=True, figsize=(700 // 72, 250 // 72))
    # ax1 = fig.add_subplot(311)
    # ax2 = fig.add_subplot(312, sharex=ax1)
    # ax3 = fig.add_subplot(313, sharex=ax1)

    ax1.plot(x * 6e0, func(x))
    ax1.plot(xs2 * 6e0, func(xs2), '.', color="C2", label="Adaptive")
    ax1.plot(x2 * 6e0, func(x2), 'x', color="C1", label="Simpson rule")
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax1.legend()

    ax2.plot(x2 * 6e0, np.zeros_like(x2), ".", color="C1")
    ax2.set_yticks([])
    plt.setp(ax2.get_xticklabels(), visible=False)

    ax3.plot(xs2 * 6e0, np.zeros_like(xs2), ".", color="C2")
    ax3.set_yticks([])
    ax3.set_xlabel("Integration variable")
    
    fig.tight_layout()

    plt.show()
