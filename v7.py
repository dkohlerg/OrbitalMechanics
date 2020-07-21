import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D
import planetary_data as pd
from OrbitPropagator import OrbitPropagator as OP
import tools as t

plt.style.use('dark_background')


# Timespan : 1 day
tspan = 3600 * 24 * 1.0
# Timestep
dt = 100.0

#central body
cb = pd.earth

if __name__ == '__main__':

    op = OP(t.tle2coes('TLE-info.txt'),tspan,dt,coes=True,deg=False)
    op.plot3d(show_plot=True)
