import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D
import planetary_data as pd
import math as m
import datetime


# plot multiple orbits
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

#convert classical orbital elements to r and v vectors
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

# convert TLE to classical orbital elements
def tle2coes(tle_filename,mu=pd.earth['mu']):
    with open(tle_filename, 'r') as f:
        lines=f.readlines()

    line0 = lines[0].strip()
    line1 = lines[1].strip().split()
    line2 = lines[2].strip().split()

    d2r = np.pi/180

    #epoch year and day
    epoch = line1[3]
    y,m,d,h = calc_epoch(epoch)

    e = float('0.' + line2[4])
    i = float(line2[2])*d2r
    aop = float(line2[5])*d2r
    raan = float(line2[3])*d2r
    Me = float(line2[6])*d2r
    mean_motion = float(line2[7]) #revolutions/day
    T = 1/(mean_motion*np.pi*2)*24*3600 #seconds
    a = (T**2*mu)**(1/3.0)

    E = ecc_anomaly([Me,e],'newton')
    ta = true_annomaly([E,e])

    return a,e,i,ta,aop,raan,[y,m,d,h]

# convert r and v vectors to classical orbital elements
def rv2coes(r,v,mu=pd.earth['mu'],deg=False,print_results=False):

    # norm of position vectors
    norm_r = np.linalg.norm(r)

    #angular momentum
    h = np.cross(r,v)
    norm_h = np.linalg.norm(h)

    # inclination
    i = m.acos(h[2]/norm_h)

    # eccenntricity vector
    e = ((np.linalg.norm(v)**2-mu/norm_r)*r-np.dot(r,v)*v)/mu

    # eccenntricity scalar
    norm_e = np.linalg.norm(e)

    # node line
    N = np.cross([0,0,1],h)
    norm_N = np.linalg.norm(N)

    # RAAN
    raan = m.acos(N[0]/norm_N)
    if N[1]<0: raan=2*np.pi-raan    #quadrant check

    # AOP
    aop = m.acos(np.dot(N,e)/(norm_N*norm_e))
    if e[2]<0: aop = 2*np.pi-aop    #quadrant check

    # True annomaly
    ta = m.acos(np.dot(e,r)/(norm_r*norm_e))
    if np.dot(r,v)<0: ta = 2*np.pi-ta    #quadrant check

    # semi major axis
    a = norm_r*(1+norm_e*m.cos(ta))/(1-norm_e**2)

    if print_results:
        print('a',a)
        print('e',norm_e)
        print('i',i*180/np.pi)
        print('RAAN',raan*180/np.pi)
        print('AOP',aop*180/np.pi)
        print('TA',ta*180/np.pi)

    if deg: return [a,norm_e,i*180/np.pi,ta*180/np.pi,aop*180/np.pi,raan*180/np.pi]
    else: return [a,norm_e,i,ta,aop,raan]



# inertial to periferical rotation matrix
def eci2perif(raan,aop,i):
 row0=[-m.sin(raan)*m.cos(i)*m.sin(aop)+m.cos(raan)*m.cos(aop),m.cos(raan)*m.cos(i)*m.sin(aop)+m.sin(raan)*m.cos(aop),m.sin(i)*m.sin(aop)]
 row1=[-m.sin(raan)*m.cos(i)*m.cos(aop)-m.cos(raan)*m.sin(aop),m.cos(raan)*m.cos(i)*m.cos(aop)-m.sin(raan)*m.sin(aop),m.sin(i)*m.cos(aop)]
 row2=[m.sin(raan)*m.sin(i),-m.cos(raan)*m.sin(i),m.cos(i)]
 return np.array([row0,row1,row2])

# calculate epoch
def calc_epoch(epoch):
       year = int('20' + epoch[:2])
       epoch = epoch[2:].split('.')
       day_of_year = int(epoch[0])-1
       hour = float('0.' + epoch[1])*24.0
       date = datetime.date(year,1,1) + datetime.timedelta(day_of_year)
       month = float(date.month)
       day = float(date.day)
       return year,month,day,hour

# calculate true annomaly
def true_annomaly(arr):
    E,e = arr
    return 2.0*m.atan(m.sqrt((1+e)/(1-e))*m.tan(E/2.0))

# calculate eccentric annomaly
def ecc_anomaly(arr,method,tol=1e-8):
    if method == 'newton':
        #newton method for iteratively finding # -*- coding: utf-8 -*-
        Me,e=arr
        if Me<m.pi/2.0: E0=Me+e/2.0
        else: E0 = Me-e
        for n in range(200):
            ratio=(E0-e*m.sin(E0)-Me)/(1-e*m.cos(E0))
            if abs(ratio)<tol:
                if n == 0: return E0
                else: return E1
            else:
                E1=E0-ratio
                E0=E1
        #did not converge
        return False
    elif method == 'tae':
        ta,e = arr
        return 2*m.atan(m.sqrt((1-e)/(1+e))*m.tan(ta/2.0))
    else:
        print('Invalid method for eccentric anomaly')
