TLE='''IRIDIUM 4[-]
1 24796U 97020E 20350.89291267 .00000062 00000-0 14911-4 0 9998
2 24796 86.3955 142.9325 0002151 94.6482 265.4960 14.35270334 236537'''

name, L1, L2 = TLE.splitlines()

import numpy as np
import matplotlib.pyplot as plt
import tools as t
import planetary_data as pd
from skyfield.api import Loader, EarthSatellite

# TLE = t.tle2coes(method='TLE',tle_filename='IRIDIUM4.txt',mu=pd.earth['mu'])
# print(TLE)


halfpi, pi, twopi = [f*np.pi for f in (0.5, 1, 2)]
degs, rads        = 180/pi, pi/180
Re                =   6378.137

load = Loader('~/Documents/fishing/SkyData')  # avoid duplicating large DE files
data = load('de421.bsp')
ts   = load.timescale(builtin=True)

minutes = np.arange(0, 200, 0.1)
times   = ts.utc(2019, 7, 23, 0, minutes)

data    = load('de421.bsp')
Earth   = data['earth']
Sun     = data['sun']
Sat     = Earth + EarthSatellite(L1, L2)

sunpos, earthpos, satpos = [thing.at(times).position.km for thing in (Sun, Earth, Sat)]


sunearth, sunsat         = earthpos-sunpos, satpos-sunpos

sunearthnorm, sunsatnorm = [vec/np.sqrt((vec**2).sum(axis=0)) for vec in (sunearth, sunsat)]

angle = np.arccos((sunearthnorm * sunsatnorm).sum(axis=0))

sunearthdistance = np.sqrt((sunearth**2).sum(axis=0))

sunsatdistance = np.sqrt((sunsat**2).sum(axis=0))

limbangle = np.arctan2(6378.137, sunearthdistance)

sunlit = []
for idx, value in enumerate(angle):
    sunlit.append(((angle[idx] > limbangle[idx]) or (sunsatdistance[idx] < sunearthdistance[idx])))

print(sunlit)
