import numpy as np
import matplotlib . pyplot as plt
from scipy . integrate import solve_ivp
import scipy.signal
import scipy.linalg
import math


def model1(t,y,Tau):
    l=1
    J=1
    d=0.5
    g=9.81
    m=9
    dxdt1=y[1]
    dxdt2=(Tau - d* y[1] - m*g*l * np.sin(y[0]))/J
    return np.array([dxdt1,dxdt2])
def model_ss(t,y,Tau):
    l=1
    J=1
    d=0.5
    g=9.81
    m=9
    u=Tau
    A = np.array([[0,1],[m*g*l/J,-d/J]])
    B = np.array([[0],[1/J]])
    x1=y[0]-np.pi
    x2=y[1]
    x_vec = np.array([[x1],[x2]])
    return np.squeeze(A @ x_vec + B * u)

def model_ricatti(t,y,Tau,K):
    l=1
    J=1
    d=0.5
    g=9.81
    m=9
    u=np.squeeze(-K @ np.array([y[0]-np.pi,y[1]]))
    dxdt1=y[1]
    dxdt2=(u - d* y[1] - m*g*l * np.sin(y[0]))/J
    return np.array([dxdt1,dxdt2])

def Zadanie2():
    l=1
    J=1
    d=0.5
    g=9.81
    m=9

    t=np.linspace(0,10,1000)
    X0 = np.array([np.pi/4,0])
    y_ss=solve_ivp(model_ss,[0,10],X0,t_eval=t,args=(0,))
    y_ivp=solve_ivp(model1,[0,10],X0,t_eval=t,args=(0,))
    plt.plot(t,y_ivp.y[0],'g-',label='x1_ivp')
    plt.plot(t,y_ivp.y[1],'y-',label='x2_ivp')
    plt.xlabel('t')
    plt.legend(loc='best')
    plt.grid()
    plt.title('Zadanie 2.1, x1(0)=[pi/4,0]')
    plt.show()
    #Przebiegi niby po częsci się pokrywają ale aproksymacja ucieka do nieskonczonosci
    #Powodem tego jest aproksymowany punkt pracy, przez który druga pochodna 
    # zmiennej stanu (przyśpieszenie) ma dodatnią część mgl/J i rośnie w nieskończoność.

    # Zadanie 2.2
    A = np.array([[0,1],[m*g*l/J,-d/J]])
    B = np.array([[0],[1/J]])
    R = np.eye(1)
    Q = np.eye(2)
    P_ricatti=scipy.linalg.solve_continuous_are(A,B,Q,R)
    K=np.linalg.inv(R) @ B.T @ P_ricatti
    X0_vec=[np.pi - 0.1, np.pi - 0.3,np.pi - 0.5,np.pi/2]
    for new_x0 in X0_vec:
        X0 = np.array([new_x0,0])
        y_ricatti=solve_ivp(model_ricatti,[0,10],X0,t_eval=t,args=(0,K))
        plt.plot(t,y_ricatti.y[0],'g-',label='x1_ricatti')
        plt.plot(t,y_ricatti.y[1],'y-',label='x2_ricatti')
        plt.xlabel('t')
        plt.legend(loc='best')
        plt.grid()
        plt.title('Zadanie 2.2, x1(0)='+str(new_x0))
        plt.show()

def model_3(t,y,u):
    J1=0.04
    J2=0.3
    m=0.5
    l=0.5
    g=9.81
    k=3
    dxdt1=y[1]
    dxdt2=-(m*g*l*np.sin(y[0]))/J1 - k/J1 * (y[0]-y[2])
    dxdt3=y[3]
    dxdt4=(k/J2) * (y[0]-y[2]) + (u/J2)
    return np.array([dxdt1,dxdt2,dxdt3,dxdt4])

def model3_ricatti(t,y,K):
    J1=0.04
    J2=0.3
    m=0.5
    l=0.5
    g=9.81
    k=3
    u=np.squeeze(-K @ y)
    dxdt1=y[1]
    dxdt2=-(m*g*l*np.sin(y[0]))/J1 - k/J1 * (y[0]-y[2])
    dxdt3=y[3]
    dxdt4=(k/J2) * (y[0]-y[2]) + (u/J2)
    return np.array([dxdt1,dxdt2,dxdt3,dxdt4])

def ricatti(t,P):
    # Reshape bo funkcja siema
    P=np.reshape(P,(4,4))
    # Parametry układu + wartości Q i R
    Q = np.eye(4)
    R = np.eye(1)
    J1=0.04
    J2=0.3
    m=0.5
    l=0.5
    g=9.81
    k=3
    # Macierze stanu zlinearyzowanego układu
    A=np.array([[0,1,0,0],[-k/J1-(m*g*l)/J1,0,k/J1,0],[0,0,0,1],[k/J2,0,-k/J2,0]])
    B=np.array([[0],[0],[0],[1/J2]])
    dpdt = -(P @ A - P@B@np.linalg.inv(R)@B.T@P + A.T@P + Q)
    return np.squeeze(np.reshape(dpdt,16))

def model3_ricatti_dyn(t,y):
    global p_sol
    J1=0.04
    J2=0.3
    m=0.5
    l=0.5
    g=9.81
    k=3
    B=np.array([[0],[0],[0],[1/J2]])
    R=np.eye(1)
    n = 1000 - int(t/t_end * 1000)
    if n < 0 or n == 0:
        P_ricatti = np.array([[p_sol.y[0,0],p_sol.y[1,0],p_sol.y[2,0],p_sol.y[3,0]],
                            [p_sol.y[4,0],p_sol.y[5,0],p_sol.y[6,0],p_sol.y[7,0]],
                            [p_sol.y[8,0],p_sol.y[9,0],p_sol.y[10,0],p_sol.y[11,0]],
                            [p_sol.y[12,0],p_sol.y[13,0],p_sol.y[14,0],p_sol.y[15,0]]]) 
    if n > 999:
        P_ricatti = np.array([[p_sol.y[0,-1],p_sol.y[1,-1],p_sol.y[2,-1],p_sol.y[3,-1]],
                            [p_sol.y[4,-1],p_sol.y[5,-1],p_sol.y[6,-1],p_sol.y[7,-1]],
                            [p_sol.y[8,-1],p_sol.y[9,-1],p_sol.y[10,-1],p_sol.y[11,-1]],
                            [p_sol.y[12,-1],p_sol.y[13,-1],p_sol.y[14,-1],p_sol.y[15,-1]]])
    else:
        P_ricatti = np.array([[p_sol.y[0,n],p_sol.y[1,n],p_sol.y[2,n],p_sol.y[3,n]],
                              [p_sol.y[4,n],p_sol.y[5,n],p_sol.y[6,n],p_sol.y[7,n]],
                            [p_sol.y[8,n],p_sol.y[9,n],p_sol.y[10,n],p_sol.y[11,n]],
                            [p_sol.y[12,n],p_sol.y[13,n],p_sol.y[14,n],p_sol.y[15,n]]])


    K = np.linalg.inv(R) @ B.T @ P_ricatti
    u=np.squeeze(-K @ y)
    dxdt1=y[1]
    dxdt2=-(m*g*l*np.sin(y[0]))/J1 - k/J1 * (y[0]-y[2])
    dxdt3=y[3]
    dxdt4=(k/J2) * (y[0]-y[2]) + (u/J2)
    return np.array([dxdt1,dxdt2,dxdt3,dxdt4])

def Zadanie3():
    J1=0.04
    J2=0.3
    m=0.5
    l=0.5
    g=9.81
    k=3
    A=np.array([[0,1,0,0],[-k/J1-(m*g*l)/J1,0,k/J1,0],[0,0,0,1],[k/J2,0,-k/J2,0]])
    B=np.array([[0],[0],[0],[1/J2]])
    C=np.array([0,0,0,1])
    D=np.array([0])
    sys= scipy.signal.StateSpace(A,B,C,D)
    t=np.linspace(0,10,1000)
    X0 = np.array([np.pi/4,0,0,0])
    y_ivp=solve_ivp(model_3,[0,10],X0,t_eval=t,args=(0,))
    t_ss,y_ss,x_ss=scipy.signal.lsim(sys,0,t,X0)
    plt.figure()
    plt.subplot(2,1,1)
    plt.title("Przebiegi z modelu nieliniowego")
    plt.plot(t,y_ivp.y[0],'g-',label='x1_ivp')
    plt.plot(t,y_ivp.y[1],'y-',label='x2_ivp')
    plt.plot(t,y_ivp.y[2],'b-',label='x3_ivp')
    plt.plot(t,y_ivp.y[3],'r-',label='x4_ivp')

    plt.xlabel('t')
    plt.legend(loc='best')
    plt.grid()
    plt.subplot(2,1,2)
    plt.title("Przebiegi z modelu zlinearyzowanego")
    plt.plot(t_ss,x_ss[:,0],'g-',label='x1_ss')
    plt.plot(t_ss,x_ss[:,1],'y-',label='x2_ss')
    plt.plot(t_ss,x_ss[:,2],'b-',label='x3_ss')
    plt.plot(t_ss,x_ss[:,3],'r-',label='x4_ss')
    plt.xlabel('t')
    plt.legend(loc='best')
    plt.grid()
    plt.show()

    # Zadanie 3.3
    
    # Parametry LQR
    R = np.eye(1)
    Q = np.eye(4)
    S = 10* np.eye(4)
    
    # P i K dla lqr z nieskonczonym horyzontem czasowym
    P_ricatti=scipy.linalg.solve_continuous_are(A,B,Q,R)
    K = np.linalg.inv(R) @ B.T @ P_ricatti
    # Obliczenie dynamicznego P
    global t_end
    t_end = 1.5
    t_ricatti = np.linspace(t_end,0,1000)
    global p_sol
    p_sol = solve_ivp(ricatti,[t_end,0],np.reshape(S,(16,)),t_eval=t_ricatti)
    
    # Nowe parametry symulacji
    t_new=np.linspace(0,t_end,1000)
    X0 = np.array([np.pi,0,np.pi/2,0])
    # Obliczenie modelu z LQR statycznym i dynamicznym 
    y_ricatti=solve_ivp(model3_ricatti,[0,t_end],X0,t_eval=t_new,args=(K,))
    y_ricatti_dyn=solve_ivp(model3_ricatti_dyn,[0,t_end],X0,t_eval=t_new)
    # Wykresy
    plt.figure()
    plt.subplot(2,1,1)
    plt.title("Przebiegi LQR")
    for i in range(4):
        plt.plot(t_new,y_ricatti.y[i],label='x'+str(i+1)+'_ricatti')
    plt.xlabel('t')
    plt.legend(loc='best')
    plt.grid()
    plt.subplot(2,1,2)
    plt.title("Przebiegi LQR z dynamicznym P")
    for i in range(4):
        plt.plot(t_new,y_ricatti_dyn.y[i],label='x'+str(i+1)+'_ricatti_dyn')
    plt.xlabel('t')
    plt.legend(loc='best')
    plt.grid()
    plt.show()


if __name__ == "__main__":
    Zadanie3()


