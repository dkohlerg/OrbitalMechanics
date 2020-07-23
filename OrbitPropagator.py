
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D
import planetary_data as pd
import tools as t


def null_perts():
    return {
            "J2":False,
            "aero":False,
            "moon_g":False,
            "Solar_grav":False
    }

class OrbitPropagator:
    def __init__(self,state0,tspan,dt,coes=False,deg=True,cb=pd.earth,perts=null_perts()):
        if coes:
            self.r0,self.v0 = t.coes2rv(state0[:6],deg,mu=cb['mu'])
        else:
            self.r0=state0[:3]
            self.v0=state0[3:]
        self.tspan = tspan
        self.dt = dt
        self.cb = cb

        # Number of steps
        self.n_steps = int(np.ceil(self.tspan/self.dt))

        # Preallocate memory
        self.ys = np.zeros((self.n_steps,6))
        self.ts = np.zeros((self.n_steps,1))

        # initial conditions
        self.ys[0] = np.concatenate((self.r0,self.v0))
        self.step = 1

        # Initiate solver
        self.solver = ode(self.diffy_q)
        self.solver.set_integrator('lsoda')
        self.solver.set_initial_value(self.ys[0],0)

        # Define perturbation dictio nary

        self.perts = perts

        self.propagateOrbit()

    def diffy_q(self,q,y):

        # unpack state
        rx,ry,rz,vx,vy,vz = y
        r = np.array([rx,ry,rz])

        # norm of the radius
        norm_r = np.linalg.norm(r)

        # two body acceleration
        a = - r * self.cb['mu']/norm_r**3


        # J2 perturbation

        if self.perts["J2"]:
            z2 = r[2]**2
            r2 = norm_r**2
            tx = r[0]/norm_r*(5*z2/r2-1)
            ty = r[1]/norm_r*(5*z2/r2-1)
            tz = r[2]/norm_r*(5*z2/r2-3)

            a_j2 = 1.5*self.cb['J2']*self.cb['mu']*self.cb['radius']**2/norm_r**4*np.array([tx,ty,tz])

            a+=a_j2

        return [vx,vy,vz,a[0],a[1],a[2]]

    def propagateOrbit(self):
        # Propagate orbit
        while self.solver.successful() and self.step<self.n_steps:
            self.solver.integrate(self.solver.t+self.dt)
            self.ts[self.step] = self.solver.t
            self.ys[self.step] = self.solver.y
            self.step+=1
        self.rs=self.ys[:,:3]
        self.vs=self.ys[:,3:]

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
