import numpy as np
import matplotlib . pyplot as plt
from scipy . integrate import solve_ivp
from scipy . integrate import odeint
import scipy.signal
import scipy.linalg
import math

def model_dyn(t,y,k1,k2):
    #parametry modelu
    C = 0.01
    L = 0.1
    R = 5
    Vin = 8
    #Zmienne modelu
    x1 = y[0]
    x2 = y[1]
    dxdt = np.array([[x1], [x2]]).reshape(2,1)
    e1 = 5-x1
    e2 = 0-x2
    dedt = np.array([[e1], [e2]]).reshape(2,1)
    # Wejścia
    uFB = k1*e1 + k2*e2
    vd, vd_prim, vd_bis = generator(t)
    uFF = 1/Vin*(C*L*vd_bis + L/R*vd_prim + vd)
    u = uFB + uFF
    #Macierze Stanu
    A = np.array([[0, 1],[-1/(C*L), -1/(C*R)]])
    B = np.array([[0],[Vin/(C*L)]])
    K = np.array([[k1,k2]])
    eprim = np.dot((A-np.dot(B,K)),dxdt) 
    # BK = np.array([[0, 0],[Vin/(C*L)*k1, Vin/(C*L)*k2]])
    xprim = A@dxdt + B*u
    return np.array([xprim[0],xprim[1],eprim[0],eprim[1]]).flatten()

def generator(t):
    vd = 3 + 2*np.sin(t)
    vd_prim = 2*np.cos(t)
    vd_bis = -2*np.sin(t)
    return np.array([vd, vd_prim, vd_bis])

def Zadanie1():
    C=0.01
    L=0.1
    R=5
    Vin=8
    omega_c = np.array([1, 5, 10])
    omega_0 = 1
    t=np.linspace(0,10,1000)
    for i in range(3):
        k1=(C*L/Vin)*(omega_c**2 - 1/(C*L))
        k2=(C*L/Vin)*(2*omega_c - 1/(C*R))
        L1 = 2*omega_0 - 1/(C*R)
        L2 = (omega_0 - 1/(R*C))**2 - 1/(C*L)
        y_sol=solve_ivp(model_dyn,[0,10],y0=[2,0,0,0],t_eval=t,args=(k1,k2,))
        plt.figure()
        plt.plot(t,y_sol.y[0],'g-',label='x1')
        plt.plot(t,y_sol.y[1],'y-',label='x2')
        plt.plot(t,y_sol.y[2],'r-',label='e1')
        plt.plot(t,y_sol.y[3],'b-',label='e2')
        plt.grid()
        plt.xlabel('t')
        plt.legend(loc='best')
        plt.title('Zadanie 1, ω = '+str(omega_c[i]))

def Zadanie2():
    C=0.01
    L=0.1
    R=5
    Vin=8
    omega_c = np.array([1, 5, 10])
    omega_0 = 1
    t=np.linspace(0,10,1000)
    for i in range(3):
        k1=(C*L/Vin)*(omega_c**2 - 1/(C*L))
        k2=(C*L/Vin)*(2*omega_c - 1/(C*R))
        L1 = 2*omega_0 - 1/(C*R)
        L2 = (omega_0 - 1/(R*C))**2 - 1/(C*L)
        y_sol=solve_ivp(model_dyn,[0,10],y0=[2,0,0,0],t_eval=t,args=(k1,k2,L1,L2,))
        plt.figure()
        plt.plot(t,y_sol.y[0],'g-',label='x1')
        plt.plot(t,y_sol.y[1],'y-',label='x2')
        plt.plot(t,y_sol.y[2],'r-',label='e1')
        plt.plot(t,y_sol.y[3],'b-',label='e2')
        plt.grid()
        plt.xlabel('t')
        plt.legend(loc='best')
        plt.title('Zadanie 1, ω = '+str(omega_c[i]))

if __name__ == "__main__":
    Zadanie1()
    plt.show()

