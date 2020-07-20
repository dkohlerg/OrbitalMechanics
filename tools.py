import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D
import planetary_data as pd

def plot_n_orbits(rs,cb,labels,show_plot=False,save_plot=False,title='Many orbits'):
    fig = plt.figure(figsize=(18,6))
    ax = fig.add_subplot(111,projection='3d')

    n=0

    for r in rs:
        # plot Trajectory
        ax.plot(r[:,0],r[:,1],r[:,2],label=labels[n],zorder=10)
        ax.plot([r[0,0]],[r[0,1]],[r[0,2]],'o',zorder=10)
        n = n +1

    # plot central body
    _u,_v= np.mgrid[0:2*np.pi:20j,0:np.pi:10j]
    _x = cb['radius'] * np.cos(_u) * np.sin(_v)
    _y = cb['radius'] * np.sin(_u) * np.sin(_v)
    _z = cb['radius'] * np.cos(_v)

    ax.plot_surface(_x,_y,_z, cmap = 'gist_earth', zorder=0)

    max_val = np.max(np.abs(r))

    ax.set_xlim([-max_val,max_val])
    ax.set_ylim([-max_val,max_val])
    ax.set_zlim([-max_val,max_val])

    ax.set_xlabel(['X (Km)'])
    ax.set_ylabel(['Y (Km)'])
    ax.set_zlabel(['Z (Km)'])

    # ax.set_aspect('equal')
    ax.set_title(title)

    plt.legend()

    if show_plot:
        plt.show()
    if save_plot:
        plt.savefig(title + 'png', dpi=300)

#convert classical orvital elements to r and v vectors
