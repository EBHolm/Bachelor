#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 26 17:38:18 2019
    bplot_tools
    TOOLS FOR PLOTTING DATA FROM C++ SCRIPTS
@author: emil
"""
from WongDataClass import WongData
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colors
import pickle
import time

# change default colors in matplotlib
import matplotlib
from cycler import cycler
matplotlib.rcParams['axes.prop_cycle'] = cycler(color='brgcymk')

from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)

from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)
        
        
def bpos(filename,g_ind):
        # create string "dataparams.txt"
    paramstr1 = filename[0:-4]
    paramstr2 = "params.txt"
    paramstr = paramstr1+paramstr2
    
    data = np.loadtxt(filename,dtype=float,delimiter=', ')
    params = np.loadtxt(paramstr,dtype=float)
    
    t = data[:,0]
    t = t*0.658     # conversion to fs
    z, u = [], []
    for i in range(1,5):
        z.append(data[:,i])
        u.append(data[:,int(i+4)])
        
        # convention
        z[i-1] = np.append(z[i-1],0)
        u[i-1] = np.append(u[i-1],0)
        
    q11, q12, q13 = data[:,9], data[:,10], data[:,11]
    q21, q22, q23 = data[:,12], data[:,13], data[:,14]
    q31, q32, q33 = data[:,15], data[:,16], data[:,17]
    q41, q42, q43 = data[:,18], data[:,19], data[:,20]
    
    q1, q2, q3, q4 = [], [], [], []
    for j in range(0,len(t)):
        q1.append(np.array([q11[j],q12[j],q13[j]]))
        q2.append(np.array([q21[j],q22[j],q23[j]]))
        q3.append(np.array([q31[j],q32[j],q33[j]]))
        q4.append(np.array([q41[j],q42[j],q43[j]]))
    q = [q1, q2, q3, q4]
    
    g = [params]
    nilarray = np.array([0, 0, 0],dtype=float)
    parametres = [g, 0, 0, 0, 0, t, 0]
    initial_conditions = [0.1, 0.5, 1, 1, nilarray, nilarray, nilarray]
    new_DATA = [z, u, q, 0, 0, 0]
    
    mDATA = []
    mDATA.append(new_DATA)
    DATA = WongData(mDATA, parametres, initial_conditions)
    # particle position/velocity plot
    plt.figure()
    plt.subplot(211)
    plt.subplot(2,1,1)
#    plt.title('Particle Positions and Velocities, g='+str(DATA.params[0][g_ind]))
    plt.ylabel('Position [$\mu$m]',fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    dataset = DATA.data
    N = 4
    times = DATA.params[5]
    
    tmax = max(times)
    
    for n in range(0,N):
        z = dataset[g_ind][0][n][1:]*0.197      # conversion to micrometres
        plt.plot(times,z,'-')
        
#    plt.ylim([-0.011, 0.011])
    plt.xlim([0, tmax])
    
    plt.subplot(2,1,2)
    plt.xlabel('Time [fs]',fontsize=20)
    plt.ylabel('Velocity [c]',fontsize=20)
    for n in range(0,N):
        plt.plot(times,dataset[g_ind][1][n][1:],'-')
    plt.ylim([-0.15, 0.15])
    plt.xlim([0, tmax])
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.tight_layout()

def bcsphere(filename,parnum):
    # create string "dataparams.txt"
    paramstr1 = filename[0:-4]
    paramstr2 = "params.txt"
    paramstr = paramstr1+paramstr2
    
    data = np.loadtxt(filename,dtype=float,delimiter=', ')
    params = np.loadtxt(paramstr,dtype=float)
    
    e=0.303
    qlen=1.04
    q1, q2, q3 = data[:,int(3*parnum+6)]/e*qlen, data[:,int(3*parnum+7)]/e*qlen, data[:,int(3*parnum+8)]/e*qlen
    print(str(int(3*parnum+6))+", "+str(int(3*parnum+7))+", "+str(int(3*parnum+8)))
    print(np.array([q1[0],q2[0],q3[0]]))
    
    # charge sphere plot, n is particle number
    r = 1
    pi = np.pi
    cos = np.cos
    sin = np.sin
    phi, theta = np.mgrid[0.0:pi:100j, 0.0:2.0*pi:100j]
    cx = r*sin(phi)*cos(theta)
    cy = r*sin(phi)*sin(theta)
    cz = r*cos(phi)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
#    ax.set_title("Charge vector, particle "+str(parnum))
    ax.plot_surface(cx, cy, cz,  rstride=1, cstride=1, color='c', alpha=0.3, linewidth=0)

    
    phi,theta=np.linspace(0,2*np.pi,100), np.pi/2
    ox = r*sin(phi)*cos(theta)
    oy = r*sin(phi)*sin(theta)
    oz = r*cos(phi)
    ax.scatter(ox,oy,oz,color='k',s=1)
#    phi,theta=np.linspace(0,2*np.pi,100), 0
#    ox = r*sin(phi)*cos(theta)
#    oy = r*sin(phi)*sin(theta)
#    oz = r*cos(phi)
#    ax.scatter(ox,oy,oz,color='k',s=1)
    phi,theta= np.pi/2,np.linspace(0,2*np.pi,100)
    ox = r*sin(phi)*cos(theta)
    oy = r*sin(phi)*sin(theta)
    oz = r*cos(phi)
    ax.scatter(ox,oy,oz,color='k',s=1)
    
    ax.scatter(q1,q2,q3,color='r',s=2,alpha=0.2)
    alength=1.05*qlen
    q_i = Arrow3D([0, q1[0]*alength], [0, q2[0]*alength], [0, q3[0]*alength], mutation_scale=20, lw=3, arrowstyle="-|>", color="b")
    q_f = Arrow3D([0, q1[-1]*alength], [0, q2[-1]*alength], [0, q3[-1]*alength], mutation_scale=20, lw=3, arrowstyle="-|>", color="m")
    ax.add_artist(q_i)
    ax.add_artist(q_f)
    ax.set_xlim([-1,1])
    ax.set_ylim([-1,1])
    ax.set_zlim([-1,1])
    ax.set_xlabel("$q_1$ [e]",fontsize=20)
    ax.set_ylabel("$q_2$ [e]",fontsize=20)
    ax.set_zlabel("$q_3$ [e]",fontsize=20)
    ax.set_aspect("equal")
    plt.tight_layout()
#    plt.axis("off")
    plt.show()
    
def bposduo(filename1,filename2):
    # create string "dataparams.txt"
    data1 = np.loadtxt(filename1,dtype=float,delimiter=', ')
    data2 = np.loadtxt(filename2,dtype=float,delimiter=', ')
    
    plt.figure()
    plt.subplot(211)
    plt.subplot(2,1,1)
    plt.ylabel('Position [$\mu$m]',fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylim([-0.011, 0.011])
    
    # assuming that they have the same time vectors!
    t = data1[:,0]*0.658
    plt.xlim([0, max(t)])
    z1, z2 = [], []
    for i in range(1,5):
        plt.subplot(2,1,1)
        z1.append(data1[:,i])
        plt.plot(t,z1[int(i-1)]*0.197,'-')
        
        plt.subplot(2,1,2)
        z2.append(data2[:,i])
        plt.plot(t,z2[int(i-1)]*0.197,'-')

    plt.subplot(2,1,2)
    plt.xlim([0, max(t)])
    plt.ylim([-0.011, 0.011])
    plt.xlabel('Time [fs]',fontsize=20)
    plt.ylabel('Position [$\mu$m]',fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.tight_layout()
    
    
    
    
# function that loads data from a Boozer script 
def boozerload(filename):
    # create string "dataparams.txt"
    paramstr1 = filename[0:-4]
    paramstr2 = "params.txt"
    paramstr = paramstr1+paramstr2
    
    data = np.loadtxt(filename,dtype=float,delimiter=', ')
    params = np.loadtxt(paramstr,dtype=float)
    
    t = data[:,0]
    z, u = [], []
    for i in range(1,5):
        z.append(data[:,i])
        u.append(data[:,int(i+4)])
        
        # convention
        z[i-1] = np.append(z[i-1],0)
        u[i-1] = np.append(u[i-1],0)
        
    q11, q12, q13 = data[:,9], data[:,10], data[:,11]
    q21, q22, q23 = data[:,12], data[:,13], data[:,14]
    q31, q32, q33 = data[:,15], data[:,16], data[:,17]
    q41, q42, q43 = data[:,18], data[:,19], data[:,20]
    
    q1, q2, q3, q4 = [], [], [], []
    for j in range(0,len(t)):
        q1.append(np.array([q11[j],q12[j],q13[j]]))
        q2.append(np.array([q21[j],q22[j],q23[j]]))
        q3.append(np.array([q31[j],q32[j],q33[j]]))
        q4.append(np.array([q41[j],q42[j],q43[j]]))
    q = [q1, q2, q3, q4]
    
    g = [params]
    nilarray = np.array([0, 0, 0],dtype=float)
    parametres = [g, 0, 0, 0, 0, t, 0]
    initial_conditions = [0.1, 0.5, 1, 1, nilarray, nilarray, nilarray]
    new_DATA = [z, u, q, 0, 0, 0]
    
    DATA = []
    DATA.append(new_DATA)
    DATASET = WongData(DATA, parametres, initial_conditions)
    DATASET.save()

# plots the length of the curve traversed by the charge vector
# naming of files should be "filename.txt","filename1.txt","filename2.txt" etc.
def bgloop(filenames,parnum):
    def len3d(x,y,z):
        return np.sum(np.sqrt(np.diff(x)**2 + np.diff(y)**2 + np.diff(z)**2))
    Qlens = []
    g_vec = []
    data_vec = []
    
    for fname in filenames:
        temp1 = fname[0:-4]
        temp2 = "params.txt"
        paramstr1 = temp1+temp2
        data_vec.append(np.loadtxt(fname,dtype=float,delimiter=', '))
        g_vec.append(np.loadtxt(paramstr1,dtype=float))
        
    for ind in range(0,len(data_vec)):
        q11, q12, q13 = data_vec[ind][:,9], data_vec[ind][:,10], data_vec[ind][:,11]
        q21, q22, q23 = data_vec[ind][:,12], data_vec[ind][:,13], data_vec[ind][:,14]
        q31, q32, q33 = data_vec[ind][:,15], data_vec[ind][:,16], data_vec[ind][:,17]
        q41, q42, q43 = data_vec[ind][:,18], data_vec[ind][:,19], data_vec[ind][:,20]
        Q = [[q11, q12, q13],[q21, q22, q23],[q31,q32,q33],[q41,q42,q43]]
        t = data_vec[ind][:,0]
        tpoints = len(t)
        
        q1, q2, q3 = [], [], []
        for j in range(0,tpoints):
            q1.append(Q[parnum][0])
            q2.append(Q[parnum][1])
            q3.append(Q[parnum][2])
        Qlens.append(round(len3d(q1,q2,q3),4))
    
    fit = np.polyfit(g_vec,Qlens,1)
    glist = np.linspace(0,max(g_vec)*1.2,100)
    plt.figure()
    plt.xlabel('g',fontsize=20)
    plt.ylabel('$L_Q$ [$e$]',fontsize=20)
    plt.plot(glist,fit[0]*glist+fit[1],'r-')
    plt.plot(g_vec,Qlens,'b.',markersize=12)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
#    plt.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
    plt.tight_layout
    
# plots the deviation between two solutions
def bdev(filename1,filename2,parnum):
    # create string "dataparams.txt"
    temp1 = filename1[0:-4]
    temp2 = "params.txt"
    paramstr1 = temp1+temp2
    data1 = np.loadtxt(filename1,dtype=float,delimiter=', ')
    params1 = np.loadtxt(paramstr1,dtype=float)
    
    temp1 = filename2[0:-4]
    temp2 = "params.txt"
    paramstr2 = temp1+temp2
    data2 = np.loadtxt(filename2,dtype=float,delimiter=', ')
    params2 = np.loadtxt(paramstr2,dtype=float)
    
    # assuming that they have the same time vectors!
    t = data1[:,0]*0.658
    z1, u1 = [], []
    z2, u2 = [], []
    for i in range(1,5):
        z1.append(data1[:,i])
        u1.append(data1[:,int(i+4)])
        z2.append(data2[:,i])
        u2.append(data2[:,int(i+4)])
        
    dev = np.abs(z1[parnum]-z2[parnum])
    
    plt.figure()
    plt.title('Deviation between two solutions, particle '+str(parnum))
    plt.xlabel('Time [fs]',fontsize=20)
    plt.ylabel('Deviation [$\mu$m]',fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.plot(t,dev*0.197,'b-')
    plt.xlim([0, max(t)])
    
# plots the pythagorean sum of the deviations
def bdevall(filename1,filename2,parnum):
    # create string "dataparams.txt"
    temp1 = filename1[0:-4]
    temp2 = "params.txt"
    paramstr1 = temp1+temp2
    data1 = np.loadtxt(filename1,dtype=float,delimiter=', ')
    params1 = np.loadtxt(paramstr1,dtype=float)
    
    temp1 = filename2[0:-4]
    temp2 = "params.txt"
    paramstr2 = temp1+temp2
    data2 = np.loadtxt(filename2,dtype=float,delimiter=', ')
    params2 = np.loadtxt(paramstr2,dtype=float)
    
    # assuming that they have the same time vectors!
    t = data1[:,0]*0.658
    z1, u1 = [], []
    z2, u2 = [], []
    for i in range(1,5):
        z1.append(data1[:,i])
        u1.append(data1[:,int(i+4)])
        z2.append(data2[:,i])
        u2.append(data2[:,int(i+4)])
        
    dev = 0
    for n in range(0,4):
        dev += np.power(z1[n]-z2[n],2)
    dev = np.sqrt(dev)/4
    
    plt.figure()
#    plt.title('Deviation between two solutions, all particles')
    plt.xlabel('Time [fs]',fontsize=20)
    plt.ylabel('Deviation [$\mu$m]',fontsize=20)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.plot(t,dev*0.197,'b-',lw=2)
    plt.xlim([0, max(t)])
    plt.tight_layout()
    
def bmultidev(refname,filenames):
    g_vec = []
    data_vec = []
    
    data_vec.append(np.loadtxt(refname,dtype=float,delimiter=', '))
    g_vec.append(0)  # start with g=0
    
    for fname in filenames:
        temp1 = fname[0:-4]
        temp2 = "params.txt"
        paramstr1 = temp1+temp2
        data_vec.append(np.loadtxt(fname,dtype=float,delimiter=', '))
        g_vec.append(np.loadtxt(paramstr1,dtype=float))
    
    devs = []
    for k in data_vec[1:]:
        dev = 0
        for n in range(0,4):
            dev += np.power(k[:,int(n+1)]-data_vec[0][:,int(n+1)],2)
        dev = np.sqrt(dev)/4
        devs.append(dev)
    
    t = data_vec[0][:,0]
    
    plt.figure()
    plt.xlabel('Time [fs]',fontsize=20)
    plt.ylabel('Deviation [$\mu$m]',fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlim([0, max(t)/2])
    for dev in devs:
        plt.plot(t,dev,'-')
    plt.tight_layout()
    plt.text(2.87,0.011,'$g=25$',fontsize=18)
    plt.text(4.455,0.004494,'$g=15$',fontsize=18)
    plt.text(7.20,0.002809,'$g=5$',fontsize=18)
    plt.ylim([0, 0.028])
    
def bcumdev(refname,filenames):
    g_vec = []
    data_vec = []
    
    data_vec.append(np.loadtxt(refname,dtype=float,delimiter=', '))
    g_vec.append(0)  # start with g=0
    
    for fname in filenames:
        temp1 = fname[0:-4]
        temp2 = "params.txt"
        paramstr1 = temp1+temp2
        data_vec.append(np.loadtxt(fname,dtype=float,delimiter=', '))
        g_vec.append(np.loadtxt(paramstr1,dtype=float))
#    cumdev = 0
    devs = []
    for k in data_vec[1:]:
        dev = 0
        for n in range(0,4):
            dev += np.power(k[:,int(n+1)]-data_vec[0][:,int(n+1)],2)
        dev = np.sqrt(dev)/4
#        cumdev += dev
        devs.append(dev)
#        devs.append(cumdev)
    
    t = data_vec[0][:,0]
    
    plt.figure()
    plt.xlabel('Time [fs]',fontsize=20)
    plt.ylabel('Cumulative Deviation [$\mu$m]',fontsize=20)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlim([0, max(t)])
    for dev in devs:
        plt.plot(t,np.cumsum(dev),'-')
    plt.tight_layout()
    plt.legend(["$\Delta \ z(0)=10^{-3} \ \mu$m",
               "$\Delta \ z(0)=10^{-4} \ \mu$m",
               "$\Delta \ z(0)=10^{-6} \ \mu$m",
               "$\Delta \ z(0)=10^{-10} \ \mu$m"],fontsize=16,frameon=False)
#    plt.text(2.87,0.011,'$g=25$',fontsize=18)
#    plt.text(4.455,0.004494,'$g=15$',fontsize=18)
#    plt.text(7.20,0.002809,'$g=5$',fontsize=18)
#    plt.ylim([0, 0.028])
    
def bstate(filename,parnum):
    data = np.loadtxt(filename,dtype=float,delimiter=', ')
    
    t = data[:,0]*0.658
    z, u = [], []
    for i in range(1,5):
        z.append(data[:,i])
        u.append(data[:,int(i+4)])
    
    plt.figure()
    plt.title('State space orbit of particle '+str(parnum))
    plt.xlabel('$z$ [eV$^{-1}$]',fontsize=20)
    plt.ylabel('$u$ [c]',fontsize=20)
    plt.plot(z[parnum],u[parnum],'b-')
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    
def bpoinduo(filename1,filename2,n,period):
    if period == 0:
#        period = 2.22635689
        period = 2.2491
#        period = 0.74
    # calculated period is 1.465 fs or 2.22635689 eV^-1
    data1 = np.loadtxt(filename1,dtype=float,delimiter=', ')
    paramstr1 = filename1[0:-4]
    paramstr2 = "params.txt"
    paramstr = paramstr1+paramstr2
    params1 = np.loadtxt(paramstr,dtype=float)
    g1 = params1
    
    data2 = np.loadtxt(filename2,dtype=float,delimiter=', ')
    paramstr1 = filename2[0:-4]
    paramstr2 = "params.txt"
    paramstr = paramstr1+paramstr2
    params2 = np.loadtxt(paramstr,dtype=float)
    g2 = params2
    
    t = data1[:,0]
    z1, u1 = [], []
    z2, u2 = [], []
    for i in range(1,5):
        z1.append(data1[:,i])
        u1.append(data1[:,int(i+4)])
        z2.append(data2[:,i])
        u2.append(data2[:,int(i+4)])
        
    # make time vector from 0 to tmax in steps of T
    tmax = max(t)
    steps = np.floor(tmax/period)
    T_times = np.linspace(0,tmax,steps)
    # find points corresponding to T_times
    poinpos1, poinvel1 = [], []
    poinpos2, poinvel2 = [], []
    t_inds = [0]
    count=1
    
    for T in T_times:
        t_ind = np.argmin(np.abs(count*period-t[int(t_inds[count-1]+10):]))
#        print(int(t_inds[count-1]+10))
        t_ind += t_inds[count-1]+10
        t_inds.append(t_ind)
        count += 1
#        print(count)
        poinpos1.append(z1[n][t_ind]*0.197)
        poinvel1.append(u1[n][t_ind])
        poinpos2.append(z2[n][t_ind]*0.197)
        poinvel2.append(u2[n][t_ind])
    
    # convert units
#    poinpos = poinpos*0.197
    
#    print(t_inds)
    plt.figure()
    plt.title('Poincaré section, particle '+str(n)+", period="+(str(period)))
    plt.xlabel('Position [$\mu$m]',fontsize=20)
    plt.ylabel('Velocity [c]',fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.plot(poinpos2,poinvel2,'b.')
    plt.plot(poinpos1,poinvel1,'r.')
    plt.tight_layout()
#    plt.figure()
#    plt.plot(t,z[n],'b.')
#    plt.plot(t[t_inds],z[n][t_inds],'r.')
#    plt.plot([0,max(t)],[0.02,0.02],'k-')
    
def bpoin(filename,n,period):
    if period == 0:
#        period = 2.22635689
        period = 2.2491
#        period = 0.74
    # calculated period is 1.465 fs or 2.22635689 eV^-1
    data = np.loadtxt(filename,dtype=float,delimiter=', ')
    paramstr1 = filename[0:-4]
    paramstr2 = "params.txt"
    paramstr = paramstr1+paramstr2
    params = np.loadtxt(paramstr,dtype=float)
    g = params
    
    
    t = data[:,0]
#    periodind = np.argmin(np.abs(t-period))
    z, u = [], []
    for i in range(1,5):
        z.append(data[:,i])
        u.append(data[:,int(i+4)])
    
    # make time vector from 0 to tmax in steps of T
    tmax = max(t)
    steps = np.floor(tmax/period)
    T_times = np.linspace(0,tmax,steps)
#    indstep = int(np.floor((T_times[1]/tmax*len(t))))    
#    print(T_times)
    # find points corresponding to T_times
    poinpos, poinvel = [], []
    t_inds = [0]
    count=1
#    print(len(t))
    
    for T in T_times:
#        print(T)
        t_ind = np.argmin(np.abs(count*period-t[int(t_inds[count-1]+10):]))
#        print(int(t_inds[count-1]+10))
        t_ind += t_inds[count-1]+10
        t_inds.append(t_ind)
        count += 1
#        print(count)
        poinpos.append(z[n][t_ind])
        poinvel.append(u[n][t_ind])
    
    # convert units
#    poinpos = poinpos*0.197
    
#    print(t_inds)
    plt.figure()
    plt.title('Poincaré section, particle '+str(n)+", period="+(str(period))+", g="+str(g))
    plt.xlabel('Position [$\mu$m]',fontsize=20)
    plt.ylabel('Velocity [c]',fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlim([-1.5, 1.5])
    plt.ylim([-0.15, 0.15])
#    plt.plot(z[n][::indstep],u[n][::indstep],'b.')
    plt.plot(poinpos,poinvel,'b.')
    plt.tight_layout
    #    plt.plot(z[n][111],u[n][111],'r.')
#    plt.plot(z[n][211],u[n][211],'r.')
#    plt.figure()
#    plt.plot(t,z[n],'b.')
#    plt.plot(t[t_inds],z[n][t_inds],'r.')
#    plt.plot([0,max(t)],[0.02,0.02],'k-')
#    plt.tight_layout()

    
def heatmap(filename,parnum,tmin,tmax):
    data = np.loadtxt(filename,dtype=float,delimiter=', ')
    Lx = 1
    N = 2001
    x = np.linspace(-Lx,Lx,N)
    
    ux = Lx
    uN = N
    xu = np.linspace(-ux,ux,uN)
    
    heatgrid = np.zeros([N, N])
    
    t = data[:,0]
    if tmax != 0:
        tmaxind = np.argmin(np.abs(t-tmax))
    else:
        tmaxind = len(t)
    if tmin != 0:
        tminind = np.argmin(np.abs(t-tmin))
    else:
        tminind = 0
    
    z, u = [], []
    for i in range(1,5):
        z.append(data[:,i])
        u.append(data[:,int(i+4)])
        
    for j in range(tminind,tmaxind):
        if parnum != 5:
            nz = np.argmin(np.abs(x-z[parnum][j]))
            nu = np.argmin(np.abs(xu-u[parnum][j]))
            heatgrid[nu,nz] += 1
        else:
            for n in range(0,4):
                nz = np.argmin(np.abs(x-z[n][j]))
                nu = np.argmin(np.abs(xu-u[n][j]))
                heatgrid[nu,nz] += 1
            
    xp= np.unravel_index(np.argmax(heatgrid, axis=None), heatgrid.shape)
    print(xp)
    print(heatgrid[xp[1],xp[0]])
    plt.figure()
    plt.title('Heat map of state space orbit, particle '+str(parnum))
    plt.xlabel('$z$ [eV$^{-1}$]',fontsize=20)
    plt.ylabel('$u$ [c]',fontsize=20)
#    plt.plot(z[parnum],u[parnum],'b-')
    plt.imshow(heatgrid, cmap='hot', interpolation='nearest')
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
#    plt.plot(xp[1],xp[0],'r.',markersize=16)
    plt.xlim([800, 1200])
    plt.ylim([800, 1200])
    
def newperiod(filename,parnum,periodmin):
    # calculated period is 1.465 fs or 2.22635689 eV^-1
    data = np.loadtxt(filename,dtype=float,delimiter=', ')
    paramstr1 = filename[0:-4]
    paramstr2 = "params.txt"
    paramstr = paramstr1+paramstr2
    params = np.loadtxt(paramstr,dtype=float)
    g = params
    
    t = data[:,0]
    x0 = np.zeros(len(t))
    z, u = [], []
    for i in range(1,5):
        z.append(data[:,i])
        u.append(data[:,int(i+4)])
    
    if periodmin == 0:
        pind = np.argmin(0-z[parnum][100:])
        print("period consists of "+str(pind)+" indices")
    else:
        pind = periodmin
    
    count = 0
    steps = int(np.floor(len(t)/pind))
    z0inds = []
    for k in range(0,steps):
        z0inds.append(np.argmin(np.abs(0-z[parnum][k*pind:])))
    print(z0inds)
    print(steps)
    print(len(t))
        
    
#    tmax = max(t)
#    steps = np.floor(tmax/period)
#    T_times = np.linspace(0,tmax,steps)
##    indstep = int(np.floor((T_times[1]/tmax*len(t))))    
##    print(T_times)
#    # find points corresponding to T_times
#    poinpos, poinvel = [], []
#    t_inds = [0]
#    count=1
##    print(len(t))
#    
#    for T in T_times:
#        t_ind = np.argmin(np.abs(count*period-t[int(t_inds[count-1]+10):]))
##        print(int(t_inds[count-1]+10))
#        t_ind += t_inds[count-1]+10
#        t_inds.append(t_ind)
#        count += 1
##        print(count)
#        poinpos.append(z[n][t_ind])
#        poinvel.append(u[n][t_ind])
#    
    
    
    
    