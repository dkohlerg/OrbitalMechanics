import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D
import planetary_data as pd
import math as m

def plot_n_orbits(rs,labels,cb=pd.earth,show_plot=False,save_plot=False,title='Many orbits'):
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
def coes2rv(coes,deg=False,mu=pd.earth['mu']):
    if deg:
        a,e,i,ta,aop,raan = coes
        i*=(np.pi/180)
        ta*=(np.pi/180)
        aop*=(np.pi/180)
        raan*=(np.pi/180)
    else:
        a,e,i,ta,aop,raan = coes

    E = ecc_anomaly([ta,e],'tae')
    r_norm = a*(1-e**2)/(1+e*m.cos(ta))

    # calculate r and v vectors in perifocal frame
    r_perif=r_norm*np.array([m.cos(ta),m.sin(ta),0])
    v_perif=m.sqrt(mu*a)/r_norm*np.array([-m.sin(E),m.cos(E)*m.sqrt(1-e**2),0])

    # rotation matrix from perifocal to ECI
    perif2eci=np.transpose(eci2perif(raan,aop,i))

    # calculate r and v vectors in inertial frames
    r = np.dot(perif2eci,r_perif)
    v = np.dot(perif2eci,v_perif)

    return r,v

# inertial to periferical rotation matrix
def eci2perif(raan,aop,i):
    row0=[-m.sin(raan)*m.cos(i)*m.sin(aop)+m.cos(raan)*m.cos(aop),m.cos(raan)*m.cos(i)*m.sin(aop)+m.sin(raan)*m.cos(aop),m.sin(i)*m.sin(aop)]
    row1=[-m.sin(raan)*m.cos(i)*m.cos(aop)-m.cos(raan)*m.sin(aop),m.cos(raan)*m.cos(i)*m.cos(aop)-m.sin(raan)*m.sin(aop),m.sin(i)*m.cos(aop)]
    row2=[m.sin(raan)*m.sin(i),-m.cos(raan)*m.sin(i),m.cos(i)]
    return np.array([row0,row1,row2])


# calculate eccentric annomaly
def ecc_anomaly(arr,method,tol=1e-8):
    if method == 'newton':
        #newton method for iteratively finding # -*- coding: utf-8 -*-
        Me,e=arr
        if Me<m.pi/2.0: E0=Me+e/2.0
        else: E0 = Me-e
        for n in range(200):
            ratio=(E0.e*m.sin(E0)-Me)/(1-e*m.cos(E0))
            if abs(ratio)<tol:
                if n == 0: return E0
                else: return E1
            else:
                E1=E0-ratio
                E0 = E1
        #did not converge
        return False
    elif method == 'tae':
        ta,e = arr
        return 2*m.atan(m.sqrt((1-e)/(1+e))*m.tan(ta/2.0))
    else:
        print('Invalid method for eccentric anomaly')
