import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D
import planetary_data as pd
from OrbitPropagator import OrbitPropagator as OP
import tools as t
from OrbitPropagator import null_perts


plt.style.use('dark_background')

# J2 effect due to the oblateness of the earth causes the mos significant perturbation
# affects raan, aop and time since perigee pasage


# Time p arams
tspan = 3600 * 24 * 10.0
dt = 100.0

#central body
cb = pd.earth

# perts

perts = null_perts()

perts['J2'] = True

if __name__ == '__main__':

    op = OP(t.tle2coes('TLE-info.txt'),tspan,dt,coes=True,deg=False,perts=perts)
    op.plot3d(show_plot=True)
