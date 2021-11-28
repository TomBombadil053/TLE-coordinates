#!/usr/bin/env python

import math as m
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import planetary_data as pd
#Personal tools

#WHOEVER READ THIS: BELOW PATH HAS TO BE CHANGE ON A PATH YOU SAVE THE FILE
#!sorry
# Love you

from sys import path
path.append('C:\\Users\\tomas\\Desktop\\WroSat\\python')

from OrbitPropagator import OrbitPropagator as OP
import tools as t


# Earth data
G_meters=6.67408e-11
G=G_meters*10**-8

earth={
        'name':'Earth',
        'mass':5.972e24,
        'mu': 5.972e24*G,
        'radius': 6378.0
}

#time parameers
tspan=3600*24*1.0
dt=10.0

#entral body
cb=pd.earth

if __name__ == '__main__':
        #a,e,i,ta,aop,raan = coes

#PLOTTING THE ORBIT
        # IRIDIUM 4[-]
        # 1 24796U 97020E 20350.89291267 .00000062 00000-0 14911 - 4 0 9998
        # 2 24796 86.3955 142.9325 0002151 94.6482 265.4960 14.35270334 236537

        # create orbit propagator instances

        print(t.tle2coes(method= 'Classic', tle_filename='IRIDIUM4.txt',mu=pd.earth['mu']))

        # WIZUALIZACJA!!!\/
        # a,e,i,ta,aop,raan = coes
        c0=[cb['radius']+7152,0.0002151, 86.3955, -1.6494059563047216, 1.6519227210860972, 2.494642733106795, 20357.95833334]
        op0=OP(c0,tspan,dt,coes=True)
        op0.propagate_orbit()
        t.plot_n_orbits([op0.rs], labels=["IRIDIUM4"], show_plot=True)


        # C0=OP(t.tle2coes('IRIDIUM4.txt'),tspan,dt,coes=True, deg=True)
        # t.plot_n_orbits([op1.rs],labels=["IRIDIUM4"], show_plot=True)
