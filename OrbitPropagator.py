
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D
import planetary_data as pd

class OrbitPropagator:
    def __init__(self,r0,v0,tspan,dt,coes=False,cb=pd.earth):
        self.r0 = r0
        self.v0 = v0
        self.tspan = tspan
        self.dt = dt
        self.cb = cb

    def diffy_q(self,q,y):

        # unpack state
        rx,ry,rz,vx,vy,vz = y
        r = np.array([rx,ry,rz])

        # norm of the radius
        norm_r = np.linalg.norm(r)

        # two body acceleration
        ax,ay,az = - r * self.cb['mu']/norm_r**3
        return [vx,vy,vz,ax,ay,az]

    def propagateOrbit(self):
            # Number of steps
            n_steps = int(np.ceil(self.tspan/self.dt))

            # Preallocate memory
            ys = np.zeros((n_steps,6))
            ts = np.zeros((n_steps,1))

            # initial conditions
            ys[0] = np.concatenate((self.r0,self.v0))
            step = 1

            # Initiate solver
            solver = ode(self.diffy_q)
            solver.set_integrator('lsoda')
            solver.set_initial_value(ys[0],0)

            # Propagate orbit
            while solver.successful() and step<n_steps:
                solver.integrate(solver.t+self.dt)
                ts[step] = solver.t
                ys[step] = solver.y
                step+=1
            self.rs=ys[:,:3]

    def plot3d(self,show_plot=False,save_plot=False,title='Orbit Representation'):
        fig = plt.figure(figsize=(18,6))
        ax = fig.add_subplot(111,projection='3d')

        # plot Trajectory
        ax.plot(self.rs[:,0],self.rs[:,1],self.rs[:,2],'w--',label='Trajectory', zorder= 10)
        ax.plot([self.rs[0,0]],[self.rs[0,1]],[self.rs[0,2]],'wo',label='Initial Position', zorder= 10)

        # plot central body
        _u,_v= np.mgrid[0:2*np.pi:20j,0:np.pi:10j]
        _x = self.cb['radius'] * np.cos(_u) * np.sin(_v)
        _y = self.cb['radius'] * np.sin(_u) * np.sin(_v)
        _z = self.cb['radius'] * np.cos(_v)

        ax.plot_surface(_x,_y,_z, cmap = 'gist_earth', zorder=0)

        max_val = np.max(np.abs(self.rs))

        ax.set_xlim([-max_val,max_val])
        ax.set_ylim([-max_val,max_val])
        ax.set_zlim([-max_val,max_val])

        ax.set_xlabel(['X (Km)'])
        ax.set_ylabel(['Y (Km)'])
        ax.set_zlabel(['Z (Km)'])

        # ax.set_aspect('equal')
        ax.set_title(title)

        plt.legend()

        if (show_plot):
            plt.show()

        if (save_plot):
            plt.savefig(title + '.png', dpi=300)
