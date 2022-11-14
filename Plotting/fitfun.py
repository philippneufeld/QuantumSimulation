import numpy as np
from scipy.special import wofz

def poly1(x, a, b):
    return a * x + b


def poly2(x, a, b, c):
    return a * x**2 + b * x**1 + c


def poly3(x, a, b, c, d):
    return a * x**3 + b * x**2 + c * x + d


