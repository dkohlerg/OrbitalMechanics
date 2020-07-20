import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D
import planetary_data as pd
from OrbitPropagator import OrbitPropagator as OP
import tools as t

plt.style.use('dark_background')


# Timespan : 1 day
tspan = 3600 * 5
# Timestep
dt = 100.0

if __name__ == '__main__':

    cb = pd.earth
    # Initial conditions
    r_mag0 = cb ['radius'] + 400.0 # Km
    v_mag0 = np.sqrt(cb ['mu']/r_mag0)

    # Initial Position and velocity
    r00 = np.array([r_mag0,0.01*r_mag0,-0.1 * r_mag0])
    v00 = np.array([0,v_mag0,0.3 * v_mag0])


    # Initial conditions
    r_mag1 = cb ['radius'] + 1500.0 # Km
    v_mag1 = np.sqrt(cb['mu']/r_mag1)*1.2

    # Initial Position and velocity
    r10 = np.array([r_mag1,0,0])
    v10 = np.array([0,v_mag1,0])

    op0 = OP(r00,v00,tspan,dt,cb)
    op1 = OP(r10,v10,tspan,dt,cb)

    op0.propagateOrbit()
    op1.propagateOrbit()
    t.plot_n_orbits([op0.rs,op1.rs],cb,labels=['ISS','random'],show_plot=True)
