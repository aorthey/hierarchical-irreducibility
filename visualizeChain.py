from numpy import dot
import numpy as np
from pylab import *
from cvxpy import *
from math import pi,cos,sin,acos,asin,atan

import matplotlib.pyplot as plt

def Rz(t):
        return np.array([[cos(t),-sin(t),0],[sin(t),cos(t),0],[0,0,1]])

def HTT(T):
        I=np.identity(3)
        return HT(I,T)

def HTR(R):
        return HT(R,np.array((0,0,0)))

def HT(R,T):
        H1 = np.vstack((R,[0,0,0]))
        H1 = np.hstack((H1,[[T[0]],[T[1]],[T[2]],[1]]))
        return H1

def robotFromConfig(r0,l0,t0,t1,s0,s1,sj1):
        colorL0 = 'k'
        colorJ1 = 'k'
        colorL1 = 'k'

        T0=np.array((-l0,0,0))
        T1=np.array((-r0,0,0))

        Tglob = np.array((x,y,0))
        Rglob = Rz(t0)
        R1=Rz(-t1)
        Hglob = HT(Rglob,Tglob)
        HT0 = HTT(T0)
        HT1 = HTT(T1)
        HR1 = HTR(R1)
        H1J = dot(Hglob,HT0)
        H12 = dot(dot(dot(Hglob,HT0),HR1),HT1)

        v = np.array((0,r0,0,1))
        v0 = dot(Hglob,v)
        v1 = dot(H12,v)

        p0 = np.array((0,0,0,1))
        p0 = dot(Hglob,p0)

        J0 = dot(H1J,p0)
        p1 = dot(H12,p0)

        c=plt.Circle((p0[0],p0[1]),s0,color=colorL0)
        ax.add_artist(c)
        c=plt.Circle((J0[0],J0[1]),sj1,color=colorJ1)
        ax.add_artist(c)
        c=plt.Circle((p1[0],p1[1]),s1,color=colorL1)
        ax.add_artist(c)

        #plt.plot([l1[0],v1[0]],[l1[1],v1[1]],'-k',linewidth=2)
        #plt.plot([l0[0],v0[0]],[l0[1],v0[1]],'-k',linewidth=2)

        plt.plot([p0[0],J0[0]],[p0[1],J0[1]],'-k',linewidth=2)
        plt.plot([p1[0],J0[0]],[p1[1],J0[1]],'-k',linewidth=2)

def robotFromConfigSimple(ax,tau,dtau,theta,length,delta):
        colorLinks = 'k'

        N = len(delta)
        angleFreeFloat = acos(np.dot(dtau,np.array((1,0,0)))/np.linalg.norm(dtau))

        Rglob = Rz(angleFreeFloat)
        Hglob = HT(Rglob,tau)

        HN = Hglob

        H = np.zeros((N,4,4))
        H[0,:,:] = HN
        for i in range(1,N):
                T0 = np.array((-length[i-1],0,0))
                R1=Rz(-theta[i-1])
                HT0 = HTT(T0)
                HR1 = HTR(R1)
                HN = dot(dot(HN,HR1),HT0)
                H[i,:,:] = HN

        p0t = np.array((0,0,0,1))
        p=np.zeros((N,4))

        for i in range(0,N):
                p[i,:] = dot(H[i,:,:],p0t)


        for i in range(0,N):
                c=plt.Circle((p[i,0],p[i,1]),delta[i],color=colorLinks)
                ax.add_artist(c)

        for i in range(0,N-1):
                plt.plot([p[i,0],p[i+1,0]],[p[i,1],p[i+1,1]],'-k',linewidth=2)

        v0 = np.array((0.2,0,0,1))
        v0 = dot(Hglob,v0)
        plt.plot([p[0,0],v0[0]],[p[0,1],v0[1]],'-k',linewidth=2)
        v0 = np.array((-2,0,0,1))
        v0 = dot(Hglob,v0)
        plt.plot([p[0,0],v0[0]],[p[0,1],v0[1]],'--k',linewidth=2)


def sweptVolumeCurvature(R,s0):
        t=linspace(-pi,pi,100)
        plt.plot((R-s0)*np.cos(t),(R-s0)*np.sin(t)+R,'-or')
        plt.plot((R+s0)*np.cos(t),(R+s0)*np.sin(t)+R,'-or')

def getNextConfig(R,li,ri,ti):
        ddsq = (R+sin(ti)*li)**2 + (cos(ti)*li)**2
        dd = math.sqrt(ddsq)
        lmbda = acos((ddsq + li*li - R*R)/(2*dd*li))
        gamma = acos((ddsq + ri*ri - R*R)/(2*dd*ri))

def visualizeChain(theta, gamma, tau, dtau, ddtau, L, delta):
        N = len(theta)

        fig=figure(1)
        ax = fig.gca()
        #ax = fig.gca(projection='3d')


        x=0.0
        y=0.0

        ###universal length and theta

        N = 6
        Lmax = 0.5
        deltaMax = 0.15
        theta = (pi/4)*np.ones((N-1,1))
        length = Lmax*np.ones((N-1,1))
        delta = deltaMax*np.ones((N,1))
        print delta

        tau = np.array((0.3,0.3,0))
        dtau = np.array((1,1,0))
        robotFromConfigSimple(ax,tau,dtau,theta,length,delta)

        ax.set_xlim((-3,3))
        ax.set_ylim((-3,3))
        ax.set_aspect('equal')
        plt.show()



if __name__ == '__main__':
        tau = np.array((0,0,0))
        dtau = np.array((1,0,0))
        ddtau = np.array((0,1,0))

        delta = np.array((1,0.5,0.5))
        L = np.array((0.5,0.5))
        visualizeChain([0,0],[0,0],tau,dtau,ddtau, L, delta)


