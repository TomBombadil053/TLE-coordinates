import numpy as np
import matplotlib.pyplot as plt
from sympy import *
from scipy.optimize import fsolve

# Variables:

ro = 6450  # kg/m^3
rw = 0.125 * 10 ** -3  # m
lw = 100 * 10 ** -3  # m
row = 8.2 * 10 ** -5  # ohm*m
cp = 465.2  # J/K*kg
e = 0.6
T0 = 273  # K
U = 80  # V
Tc = 2.725  # K

# Irradiance
Qo = 1366  # W/m^2

# Stefan constant
sigma = 5.67 * 10 ** -8  # W/m^2/K^4

# Geometry
Ac = np.pi * rw ** 2  # m^2 - Area of a circle

Acs = 2 * rw * lw  # m^2 - Cross-section area

At = ((2 * np.pi * rw) * lw) + ((np.pi * rw ** 2) * 2)  # m^2 - Total area

V = np.pi * rw ** 2 * lw  # m^3 - Volume

# Resistivity

R = row * lw / np.pi / rw ** 2
# print("Resistance:", R)

print("Initial assumptions\n"
      "\tDiameter: %1.5f m\n" % (2 * rw),
      "\tLength: %1.1f m\n" % lw,
      "\tDensity: %1.1f kg/m^3\n" % ro,
      "\tResistivity: %1.6f ohm*m\n" % row,
      "\tEmissivity: %1.1f\n" % e,
      "\tInitial temperature: %1.1f K\n" % T0,
      "\tMax voltage: %1.1f V\n" % U,
      "\tSpecific heat: %1.1f J/K*kg\n" % cp
      )

print("Enter final temperature (in kelvin):")
Tf = float(input())

def get_bool(prompt):
    while True:
        try:
           return {"true":True,"false":False}[input(prompt).lower()]
        except KeyError:
           print("Invalid input please enter True or False!")

respuesta = get_bool("Is the satellite in the shadow? (True/False):")
if (respuesta == True):
    QI = 0
    Pc= ((sigma * e * At * (Tc ** 4 - Tf ** 4)) + (U ** 2 / R)+(QI*Acs))/((Tf - T0) * ro * V * cp)
else:
    QI = 1366
    Pc= ((sigma * e * At * (Tc ** 4 - Tf ** 4)) + (U ** 2 / R)+(QI*Acs))/((Tf - T0) * ro * V * cp)
print(QI)
print(Pc)
# print("Enter mean motion (in orbits/day)") #parameter used to see if niuniuś is in the sunlight
# Mm=float(input())

# Heat gain and loss in one second

# QJ=U**2/R #Joule heat
# print("Joule heat:", QJ)
#
# QR=sigma*e*At*(Tc**4-Tf**4) #Radiation
# print("Radiation heat:", QR)
#
# QI=Qo*Acs #Irradiance
# print("Irradiance heat:", QI)

# t=((Tf-T0)*ro*V*cp)/(QR)
# print("Wire was heated in %1.8f s" %t)

# Numerical solution                 ##  Actually above code would be enough in simplest case but

def numerical_solution(p):
    QJ, QR, t, QI = p
    return((U ** 2 / R) - QJ, (sigma * e * At * (Tc ** 4 - Tf ** 4)) - QR, QI, ((((Tf - T0) * ro * V * cp)*t - (QR + QJ))  ))
#     else:
#         QI=1366
#         return ((U ** 2 / R) - QJ, (sigma * e * At * (Tc ** 4 - Tf ** 4)) - QR, Qo*Acs - QI, ((((Tf - T0) * ro * V * cp)*t - (QR + QJ + QI))  ))
#
QJ, QR, t, QI = fsolve(numerical_solution, (40, -0.1, 1.0, 0.0))

print("\tJoule heat:% 1.5f\n\tStefan loss: %1.5f \n\tTime: %1.5f \n\tSun heat: %1.1f \n\tJointly: %1.5f" % (QJ, QR, t,QI, QJ+QR+QI))



# if S == 0:
#     QI = 0
#     t=((Tf - T0) * ro * V * cp) /((sigma * e * At * (Tc ** 4 - Tf ** 4)) + (U ** 2 / R))
# else:
#     QI = 1366
#     return ((U ** 2 / R) - QJ, (sigma * e * At * (Tc ** 4 - Tf ** 4)) - QR, Qo * Acs - QI,
#             ((((Tf - T0) * ro * V * cp) * t - (QR + QJ + QI))))
#
# QJ, QR, t, QI = fsolve(numerical_solution, (40, -0.1, 1.0, 0.0))
#
# print("\tJoule heat:% 1.5f\n\tStefan loss: %1.5f \n\tTime: %1.5f \n\tSun heat: %1.1f \n\tJointly: %1.5f" % (
# QJ, QR, t, QI, QJ + QR + QI))
