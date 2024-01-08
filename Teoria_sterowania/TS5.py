import numpy as np
import matplotlib . pyplot as plt
from scipy . integrate import solve_ivp
import scipy.signal
import scipy.linalg
import math

def model(t,y,u):
    C=0.01
    L=0.1
    R=5
    Vin=8
    x1=y[0]
    x2=y[1]
    dxdt1=x1
    dxdt2=-1/(C*R) * x2 - 1/(C*L) * x1 + 1/(C*L) * Vin * u
    return np.array([dxdt1,dxdt2])

def Zadanie2():
    C=0.01
    L=0.1
    R=5
    Vin=8
    A = np.array([[0,1],[-1/(C*L),-1/(C*R)]])
    B = np.array([[0],[Vin/(C*L)]])
    R = np.array([[B,A@B]])
    rank = np.linalg.matrix_rank(R)
    print(rank)
    if rank == 2:
        print("System jest sterowalny")
    else:
        print("System nie jest sterowalny")

def model_dyn(t,y,k1,k2,ud):
    #Zmienna e
    C=0.01
    L=0.1
    R=5
    Vin=8
    x1=y[0]
    x2=y[1]
    u = y@[k1,k2] + ud
    print(u)
    dxdt1=x1
    dxdt2=-1/(C*R) * x2 - 1/(C*L) * x1 + 1/(C*L) * Vin * u
    return np.array([dxdt1,dxdt2])

def Zadanie3():
    C=0.01
    L=0.1
    R=5
    Vin=8
    Vd=5
    omega = np.array([-1,1,5,10])
    ud = Vd/Vin
    t=np.linspace(0,20,1000)
    for i in range(4):
        omega_c=omega[i]
        print('siema')
        print(omega_c)
        k1=(omega_c**2 - 1000)/8000
        k2=-(omega_c + 10)/4000
        result = solve_ivp(model_dyn,[0,20],[2,0],t_eval=t,args=(k1,k2,ud,))
        plt.figure()
        plt.plot(result.t, result.y[0], label='x1(t)')
        plt.plot(result.t, result.y[1], label='x2(t)')
        plt.xlabel('t')
        plt.legend(loc='best')
        plt.grid()
        plt.title('Zadanie 3, omega_c = ' + str(omega_c))
        plt.savefig('Zadanie3_'+str(i)+'.jpg')
        plt.show()
        #Stabilizacja na 5 volt

Zadanie3()
    