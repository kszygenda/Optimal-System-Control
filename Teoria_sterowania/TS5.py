import numpy as np
import matplotlib . pyplot as plt
from scipy . integrate import solve_ivp
from scipy . integrate import odeint
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
    mnozenie = A@B
    R = np.array[[B[0], mnozenie[0]],[B[1], mnozenie[1]]]
    rank = np.linalg.matrix_rank(R)
    print(rank)
    if rank == 2:
        print("System jest sterowalny")
    else:
        print("System nie jest sterowalny")

def model_dyn(y,t,k1,k2):
    #parametry modelu
    C = 0.01
    L = 0.1
    R = 5
    Vin = 8
    #Zmienne modelu
    x1 = y[0]
    x2 = y[1]
    dxdt = np.array([[x1], [x2]])
    e1 = 5-x1
    e2 = 0-x2
    dedt = np.array([[e1], [e2]])
    # Wejścia
    uFB = k1*e1 + k2*e2
    uD = 5/Vin
    u = uFB + uD
    #Macierze Stanu
    A = np.array([[0, 1],[-1/(C*L), -1/(C*R)]])
    B = np.array([[0],[Vin/(C*L)]])
    K = np.array([[k1,k2]])
    eprim = (A-B@K)@dedt
    # BK = np.array([[0, 0],[Vin/(C*L)*k1, Vin/(C*L)*k2]])
    xprim = A@dxdt + B*u
    return np.array([xprim[0],xprim[1],eprim[0],eprim[1]]).flatten()

def Zadanie3():
    C=0.01
    L=0.1
    R=5
    Vin=8
    ud = 5 / Vin
    omega = np.array([-1,1,5,10])
    # k1_vec=[-0.124875,-0.,-0.1125,-0.124875]
    # k2_vec=[-0.00225,-0.00125,0,-0.00275]
    t=np.linspace(0,10,1000)
    for i in range(4):
        omega_c=omega[i]
        k1=(C*L/Vin)*(omega_c**2 - 1/(C*L))
        k2=(C*L/Vin)*(2*omega_c - 1/(C*R))
        print(f'k1=',k1)
        print(f'k2=',k2)
        result = odeint(model_dyn,y0=[2,0,3,0],t=t,args=(k1,k2))
        plt.figure()
        plt.plot(t, result[:,0], label='x1(t)')
        plt.plot(t, result[:,1], label='x2(t)')
        plt.plot(t,result[:,2], label = 'e1(t)')
        plt.plot(t, result[:, 3], label='e2(t)')
        plt.plot()
        plt.xlabel('t')
        plt.legend(loc='best')
        plt.grid()
        plt.title('Zadanie 3, omega_c = ' + str(omega_c))
        plt.savefig('Zadanie3_'+str(i)+'.jpg')
        plt.show()
        #Stabilizacja na 5 volt

def model_dyn_Zpomiarowe(y,t,k1,k2):
    #parametry modelu
    C = 0.01
    L = 0.1
    R = 5
    Vin = 8
    #Zmienne modelu
    x1 = y[0]
    x2 = y[1]
    dxdt = np.array([[x1], [x2]])
    e1 = y[2]
    e2 = y[3]
    dedt = np.array([[e1], [e2]])
    # Wejścia
    uFB = k1*e1 + k2*e2 + k1*np.sin(4*np.pi*t)**2 + k2*np.cos(4*np.pi*t)**2
    uD = 5/Vin
    u = uFB + uD
    #Macierze Stanu
    A = np.array([[0, 1],[-1/(C*L), -1/(C*R)]])
    B = np.array([[0],[Vin/(C*L)]])
    K = np.array([[k1,k2]])
    eprim = (A-B@K)@dedt
    # BK = np.array([[0, 0],[Vin/(C*L)*k1, Vin/(C*L)*k2]])
    xprim = A@dxdt + B*u
    return np.array([xprim[0],xprim[1],eprim[0],eprim[1]]).flatten()

def Zadanie41():
    C=0.01
    L=0.1
    R=5
    Vin=8
    ud = 5 / Vin
    omega = np.array([1,2,5])
    # k1_vec=[-0.124875,-0.,-0.1125,-0.124875]
    # k2_vec=[-0.00225,-0.00125,0,-0.00275]
    t=np.linspace(0,10,1000)
    for i in range(3):
        omega_c=omega[i]
        k1=(C*L/Vin)*(omega_c**2 - 1/(C*L))
        k2=(C*L/Vin)*(2*omega_c - 1/(C*R))
        print(f'k1=',k1)
        print(f'k2=',k2)
        result = odeint(model_dyn_Zpomiarowe,y0=[0,0,0,0],t=t,args=(k1,k2))
        plt.figure()
        plt.plot(t, result[:,0], label='x1(t)')
        plt.plot(t, result[:,1], label='x2(t)')
        plt.plot(t, 5 - result[:,0], label = 'e1(t)')
        plt.plot(t, 0 - result[:, 1], label='e2(t)')
        plt.plot()
        plt.xlabel('t')
        plt.legend(loc='best')
        plt.grid()
        plt.title('Zadanie 4.1, omega_c = ' + str(omega_c))
        plt.savefig('Zadanie41_'+str(i)+'.jpg')
        plt.show()
        #Stabilizacja na 5 volt

def model_dyn_Zprocesu(y,t,k1,k2):
    #parametry modelu
    C = 0.01
    L = 0.1
    R = 5
    Vin = 8
    #Zmienne modelu
    x1 = y[0]
    x2 = y[1]
    dxdt = np.array([[x1], [x2]])
    e1 = y[2]
    e2 = y[3]
    dedt = np.array([[e1], [e2]])
    # Wejścia
    uFB = k1*e1 + k2*e2
    uD = 0
    u = uFB + uD + 0.5*np.sin(0.2*np.pi*t)**2
    #Macierze Stanu
    A = np.array([[0, 1],[-1/(C*L), -1/(C*R)]])
    B = np.array([[0],[Vin/(C*L)]])
    K = np.array([[k1,k2]])
    eprim = (A-B@K)@dedt
    # BK = np.array([[0, 0],[Vin/(C*L)*k1, Vin/(C*L)*k2]])
    xprim = A@dxdt + B*u
    return np.array([xprim[0],xprim[1],eprim[0],eprim[1]]).flatten()

def Zadanie42():
    C=0.01
    L=0.1
    R=5
    Vin=8
    ud = 5 / Vin
    omega = np.array([1,2,5])
    # k1_vec=[-0.124875,-0.,-0.1125,-0.124875]
    # k2_vec=[-0.00225,-0.00125,0,-0.00275]
    t=np.linspace(0,10,1000)
    for i in range(3):
        omega_c=omega[i]
        k1=(C*L/Vin)*(omega_c**2 - 1/(C*L))
        k2=(C*L/Vin)*(2*omega_c - 1/(C*R))
        print(f'k1=',k1)
        print(f'k2=',k2)
        result = odeint(model_dyn_Zprocesu,y0=[0,0,0,0],t=t,args=(k1,k2))
        plt.figure()
        plt.plot(t, result[:,0], label='x1(t)')
        plt.plot(t, result[:,1], label='x2(t)')
        plt.plot(t,0 - result[:,0], label = 'e1(t)')
        plt.plot(t, 0 - result[:, 1], label='e2(t)')
        plt.plot()
        plt.xlabel('t')
        plt.legend(loc='best')
        plt.grid()
        plt.title('Zadanie 4.2, omega_c = ' + str(omega_c))
        plt.savefig('Zadanie42_'+str(i)+'.jpg')
        plt.show()
        #Stabilizacja na 5 volt

if __name__ == "__main__":
    Zadanie3()
    Zadanie41()
    Zadanie42()
    