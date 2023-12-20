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
    dxdt2=1/J*u - d/J*y[1] - ((m*g*l)/J)*np.sin(y[0])
    return np.array([dxdt1,dxdt2])

def Zadanie51():
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
        y=solve_ivp(model51,[0,10],y0,t_eval=t,args=(u,))
        u_new=u*u_ss
        tss,yss,xss=scipy.signal.lsim(ss,u_new,t,y0)
        plt.figure()
        plt.title('u='+str(u))
        plt.plot(y.t,y.y[0],'b-',label='x1')
        plt.plot(y.t,y.y[1],'r-',label='x2')
        plt.plot(tss,xss[:,0],'g-',label='x1_ss')
        plt.plot(tss,xss[:,1],'y-',label='x2_ss')
        plt.legend(loc='best')
        plt.xlabel('t')
        plt.grid()
        plt.show()


def Zadanie52():
        # Linearyzacja innego punktu pracy
    # Dla tych punktów pracy tworzymy nowe macierze stanu.
    # pochodne będą te same ale wartości początkowe są inne
    t=np.linspace(0,10,1000)
    u_ss=np.ones_like(t)
    J=1
    g=10
    d=0.5
    m=9
    l=1
    u0=45*np.sqrt(2)
    y0=np.array([np.pi/4,0])
    # Zmiana w a polega na nowej wartości wewnątrz cosinusa
    A=np.array([[0,1],[-m*g*l/J * np.cos(y0[0]), -d/J]])
    # W B nie ma zmiany bo pochodna po u dalej ta sama powinna być, innej wartości początkowej
    B=np.array([[0],[1/J]])
    # C i D to wyjscie wiec nie ma potrzebyyn zmiany
    C=np.array([[1,0],[0,1]])
    D=np.array([[0],[0]])
    ss_new=scipy.signal.StateSpace(A,B,C,D)
    u_vec = [45*np.sqrt(2), 45*np.sqrt(2) + 2,45*np.sqrt(2) + 10,45*np.sqrt(2) + 30]
    for u in u_vec:
        print(f'U={u}')
        u_solve=u
        y0=np.array([np.pi/4,0])
        y=solve_ivp(model51,[0,10],y0,t_eval=t,args=(u_solve,))
        u_new=u*u_ss - u0
        tss,yss,xss=scipy.signal.lsim(ss_new,u_new,t,y0)
        plt.figure()
        plt.title('u='+str(u))
        plt.plot(y.t,y.y[0],'b-',label='x1')
        plt.plot(y.t,y.y[1],'r-',label='x2')
        plt.plot(tss,xss[:,0] - np.pi/4,'g-',label='x1_ss')
        plt.plot(tss,xss[:,1],'y-',label='x2_ss')
        plt.legend(loc='best')
        plt.xlabel('t')
        plt.grid()
        plt.show()

def SDC_Model(t,y,u):
    # Parametry statyczne
    J=1
    g=10
    d=0.5
    m=9
    l=1
    # Zmienne i macierze SDC 
    x1=y[0]
    x2=y[1]
    A = np.array([[0,1],[-m*g*l/J*(np.sin(x1)/x1),-d/J]])
    B = np.array([[0],[1/J]])
    x_vec = np.array([[x1],[x2]])
    return np.squeeze(A @ x_vec + B * u)

def Zadanie53():
    t=np.linspace(0,10,1000)
    x1_vec=[np.pi/4,0]
    for x1 in x1_vec:
        y0=np.array([x1,0])
        y_sol=solve_ivp(model51,[0,10],y0,t_eval=t,args=(0,))
        y_sdc=solve_ivp(SDC_Model,[0,10],y0,t_eval=t,args=(0,))
        plt.figure()
        plt.title('x1(0)='+str(x1))
        plt.plot(y_sol.t,y_sol.y[0],'b-',label='x1_ivp')
        plt.plot(y_sol.t,y_sol.y[1],'r-',label='x2_ivp')
        plt.plot(y_sdc.t,y_sdc.y[0],'g-',label='x1_sdc')
        plt.plot(y_sdc.t,y_sdc.y[1],'y-',label='x2_sdc')
        plt.legend(loc='best')
        plt.xlabel('t')
        plt.grid()
        plt.show()

if __name__ == '__main__':
    #Zadanie2()
    Zadanie53()