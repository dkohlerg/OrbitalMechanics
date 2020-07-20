import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D
import planetary_data as pd
from OrbitPropagator import OrbitPropagator as OP
import tools as t

plt.style.use('dark_background')


# Timespan : 1 day
tspan = 3600 * 24
# Timestep
dt = 10.0

#central body
cb = pd.earth

if __name__ == '__main__':

    #a,e,i,ta,aop,raan=coes

    #ISS
    c0 =[cb['radius']+414.0,0.0006189,51.6393,0.0,234.1955,105.6372]

    #GEO
    c1=[cb['radius']+35800.0,0.0,0.0,0.0,0.0,0.0]


    op0 = OP(c0,tspan,dt,coes=True)
    op1 = OP(c1,tspan,dt,coes=True)

    op0.propagateOrbit()
    op1.propagateOrbit()
    t.plot_n_orbits([op0.rs,op1.rs],labels=['ISS','GEO'],show_plot=True)
