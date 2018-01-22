# Water content to permittivity

import math


eps_w = 80.1
eps_s = 2.23
eps_a = 1


def CRIM(theta, poro):
    eps = (math.sqrt(eps_w)*theta+(1-poro)*math.sqrt(eps_s)+(poro-theta))**2
    return eps


def Topp(theta):
    eps = 3.03 + 9.3*theta + 146*math.pow(theta, 2)-76.7*math.pow(theta, 3)
    return eps
