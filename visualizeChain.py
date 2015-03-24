from numpy import dot
import numpy as np
from pylab import *
from cvxpy import *
from math import pi,cos,sin,acos,asin,atan
from mpl_toolkits.mplot3d import Axes3D

import matplotlib.pyplot as plt

colorLinks = 'k'

def Rz(t):
        return np.array([[cos(t),-sin(t),0],[sin(t),cos(t),0],[0,0,1]])
def Ry(t):
        return np.array([[cos(t),0,sin(t)],[0,1,0],[-sin(t),0,cos(t)]])

def HTT(T):
        I=np.identity(3)
        return HT(I,T)

def HTR(R):
        return HT(R,np.array((0,0,0)))

def HT(R,T):
        H1 = np.vstack((R,[0,0,0]))
        H1 = np.hstack((H1,[[T[0]],[T[1]],[T[2]],[1]]))
        return H1


def drawSphere(ax,x,y,z,R):
        u, v = np.mgrid[0:2*np.pi:10j, 0:np.pi:10j]
        xx=np.cos(u)*np.sin(v)
        yy=np.sin(u)*np.sin(v)
        zz=np.cos(v)
        xx = R*xx + x
        yy = R*yy + y
        zz = R*zz + z
        ax.plot_wireframe(xx, yy, zz, color=colorLinks)

def robotFromConfigSimple(ax,tau,dtau,ddtau,theta,gamma,length,delta):

        N = len(delta)
        angleFreeFloat = acos(np.dot(dtau,np.array((1,0,0)))/np.linalg.norm(dtau))

        Rglob = Rz(angleFreeFloat)
        Hglob = HT(Rglob,tau)

        HN = Hglob

        H = np.zeros((N,4,4))
        H[0,:,:] = HN
        for i in range(1,N):
                T0 = np.array((-length[i-1],0,0))
                Rtheta=Rz(-theta[i-1])
                Rgamma=Ry(-gamma[i-1])
                HT0 = HTT(T0)
                HRZ = HTR(Rtheta)
                HRY = HTR(Rgamma)
                HN = dot(dot(dot(HN,HRZ),HRY),HT0)
                H[i,:,:] = HN

        p0t = np.array((0,0,0,1))
        p=np.zeros((N,4))

        for i in range(0,N):
                p[i,:] = dot(H[i,:,:],p0t)


        for i in range(0,N):
                c=plt.Circle((p[i,0],p[i,1]),delta[i],color=colorLinks)
                ax.add_artist(c)

        for i in range(0,N-1):
                plt.plot([p[i,0],p[i+1,0]],[p[i,1],p[i+1,1]],[p[i,2],p[i+1,2]],'-k',linewidth=2)

        print dtau
        v0 = np.hstack([0.5*dtau,1])
        v0 = dot(Hglob,v0)
        plt.plot([p[0,0],v0[0]],[p[0,1],v0[1]],'-k',linewidth=2)
        v0 = np.hstack([-2*dtau,1])
        v0 = dot(Hglob,v0)
        plt.plot([p[0,0],v0[0]],[p[0,1],v0[1]],'--k',linewidth=2)

        ax.scatter( p[:,0],p[:,1],p[:,2],s=100)
        
        for i in range(0,N):
                drawSphere(ax, p[i,0],p[i,1],p[i,2],delta[i])


def visualizeChain(theta, gamma, tau, dtau, ddtau, L, delta):
        N = len(theta)

        fig=figure(1)
        ax = fig.gca(projection='3d')

        x=0.0
        y=0.0

        ###universal length and theta

        N = 6
        Lmax = 0.5
        deltaMax = 0.15
        theta = (pi/4)*np.ones((N-1,1)).flatten()
        gamma = (pi/8)*np.ones((N-1,1)).flatten()
        length = Lmax*np.ones((N-1,1)).flatten()
        delta = deltaMax*np.ones((N,1)).flatten()

        tau = np.array((0.3,0.3,0))
        dtau = np.array((1,0,0))
        ddtau = np.array((0,1,0))

        robotFromConfigSimple(ax,tau,dtau,ddtau,theta,gamma,length,delta)

        lim = 2
        ax.set_xlim((-lim,lim))
        ax.set_ylim((-lim,lim))
        ax.set_zlim((-lim,lim))
        ax.set_aspect('equal')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        plt.show()



if __name__ == '__main__':
        tau = np.array((0,0,0))
        dtau = np.array((1,0,0))
        ddtau = np.array((0,1,0))

        delta = np.array((1,0.5,0.5))
        L = np.array((0.5,0.5))
        visualizeChain([0,0],[0,0],tau,dtau,ddtau, L, delta)


