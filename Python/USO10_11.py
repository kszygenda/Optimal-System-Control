import numpy as np
import matplotlib . pyplot as plt
from scipy . integrate import solve_ivp
import scipy.signal
import math

def model1(t,y):
    u=1
    dzdt1 = y[0]*np.log(y[1])
    dzdt2 = -y[1]*np.log(y[0]) + y[1]*u
    return np.array([dzdt1,dzdt2])

def model12(t,y):
    u=1
    dxdt1=y[1]
    dxdt2=-y[0]+u
    return np.array([dxdt1,dxdt2])

def Zadanie2():
    t = np.linspace(0, 10, 1000)
    y0 = np.array([1, 1])
    y = solve_ivp(model1, [0, 10], y0, t_eval=t)
    y2= solve_ivp(model12,[0,10],y0=np.array([1,1]),t_eval=t)
    plt.plot(y.t, y.y[0], 'b-', label='z1')
    plt.plot(y.t, y.y[1], 'r-', label='z2')
    plt.plot(y2.t,y2.y[0],'g-',label='x1')
    plt.plot(y2.t,y2.y[1],'y-',label='x2')
    plt.xlabel('t')
    plt.legend(loc='best')
    plt.grid()
    plt.show()

def model51(t,y,u):
    J=1
    g=10
    d=0.5
    m=9
    l=1
    dxdt1=y[1]
    dxdt2=1/J*u - d/J*y[1] - m*g*l/J*np.sin(y[0])
    return np.array([dxdt1,dxdt2])
def Zadanie5():
    J=1
    g=10
    d=0.5
    m=9
    l=1
    A=np.array([[0,1],[-m*g*l/J, -d/J]])
    B=np.array([[0],[1/J]])
    C=np.array([[1,0],[0,1]])
    D=np.array([[0],[0]])
    ss=scipy.signal.StateSpace(A,B,C,D)
    t=np.linspace(0,10,1000)
    u_ss=np.ones_like(t)
    y0=np.array([0,0])  
    u_vec = [0,5,20,45*np.sqrt(2),70]
    for u in u_vec:
        y=solve_ivp(model51,[0,10],y0,t_eval=t,args=(1,))
        u_new=u*u_ss
        y2=scipy.signal.lsim(ss,u_new,t,y0)
        plt.figure
        
        plt.plot(y.t,y.y[0],'b-',label='x1')
        plt.plot(y.t,y.y[1],'r-',label='x2')
        plt.plot(y2.t,y2.y[0],'g-',label='x1_ss')
        plt.plot(y2.t,y2.y[1],'y-',label='x2_ss')
        plt.legend(loc='best')
        plt.xlabel('t')
        plt.grid()
        plt.show()

if __name__ == '__main__':
    Zadanie2()
    Zadanie5()