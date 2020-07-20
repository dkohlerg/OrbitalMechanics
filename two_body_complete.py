import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D
import planetary_data as pd
from OrbitPropagator import OrbitPropagator as OP

plt.style.use('dark_background')

cb = pd.earth

if __name__ == '__main__':
    # Initial conditions
    r_mag = cb['radius'] + 1500.0 # Km
    v_mag = np.sqrt(cb['mu']/r_mag)

    # Initial Position and velocity
    r0 = np.array([r_mag,0.01*r_mag,-0.1 * r_mag])
    v0 = np.array([0,v_mag,0.3 * v_mag])

    # Timespan : 1 day
    tspan = 3600 * 5

    # Timestep
    dt = 100.0


    op = OP(r0,v0,tspan,dt,cb)

    op.propagateOrbit()
    op.plot3d(show_plot=True)
