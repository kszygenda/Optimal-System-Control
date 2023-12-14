import numpy as np
import matplotlib . pyplot as plt
from scipy . integrate import solve_ivp
import scipy.signal
import math

def Zadanie1():
    # ------------------------------- Uklad 1 ----------------------------------
    #Dane układu
    C1=1
    C2=0.5
    #Macierze stanu 1
    A1=np.array([[-1/(2*C1),0],[0,-1/(4*C2)]])
    B1=np.array([[1/(2*C1)],[1/(4*C2)]])
    C1=np.array([0, 1])
    # Kalman1=[B1,np.dot(A1,B1)]
    # print(np.dot(A1,B1))
    # print("Macierz kalmana powinna mieć rank(K)=2")
    # print(np.linalg.matrix_rank(Kalman1))
    #Statespace 1
    SS1=scipy.signal.StateSpace(A1,B1,C1)
    Kalman1=np.transpose(np.squeeze([B1, A1 @ B1],axis=2))
    print('n=2, rank(k)=')
    print(np.linalg.matrix_rank(Kalman1))
    # ------------------------------- Uklad 2 ----------------------------------
    C1=1
    C2=2
    C3=3
    #Macierze stanu 2
    A2=np.array([[-1/C1,0,0],[0,-1/C2,0],[0,0,-1/C3]])
    B2=np.array([[1/C1],[1/C2],[1/C3]])
    C2=np.array([0, 0, 1])
    # Kalman2=[B2,np.dot(A2,B2),np.dot(np.dot(A2,A2),B2)]
    # print("Macierz kalmana powinna mieć rank(K)=3")
    # print(np.linalg.matrix_rank(Kalman2))
    SS2=scipy.signal.StateSpace(A2,B2,C2)
    Kalman2=np.transpose(np.squeeze([B2, A2 @ B2, A2@A2@B2],axis=2))
    print('n=3, rank(k)=')
    print(np.linalg.matrix_rank(Kalman2))
    # ------------------------------- Uklad 3 ----------------------------------
    L1=0.1
    L2=0.1
    C1=0.1
    C2=0.1
    R=1
    #Macierze stanu 3
    A3=np.array([[0,1/L1,0,0],[-1/(R*C1),-1/(R*C1),0,-1/(R*C1)],[0,0,0,1/L2],[0,-1/(R*C2),-1/(R*C2),-1/(R*C2)]])
    B3=np.array([[0],[-1/(R*C1)],[0],[-1/(R*C2)]])
    C3=np.array([0, 0, 0, 1])
    # Kalman3=[B3,np.dot(A3,B3),np.dot(np.dot(A3,A3,),B3),np.dot(np.dot(np.dot(A3,A3),A3),B3)]
    # print("Macierz kalmana powinna mieć rank(K)=4")
    # print(np.linalg.matrix_rank(Kalman3))
    SS3=scipy.signal.StateSpace(A3,B3,C3)
    Kalman3=np.transpose(np.squeeze([B3, A3 @ B3, A3@A3@B3,A3@A3@A3@B3],axis=2))
    print('n=4, rank(k)=')
    print(np.linalg.matrix_rank(Kalman3))
    # ------------------------------- Uklad 4 ----------------------------------
    L1=0.5
    L2=1
    C1=2
    A4=np.array([[-2/L1,0,-1/L1],[0,0,1/L2],[1/C1,-1/C1,-1/C1]])
    B4=np.array([[1/L1],[0],[0]])
    C4=np.eye(3)
    # Kalman4=[B4,np.dot(A4,B4),np.dot(np.dot(A4,A4),B4)]
    # print("Macierz kalmana powinna mieć rank(K)=3")
    # print(np.linalg.matrix_rank(Kalman4))
    SS4=scipy.signal.StateSpace(A4,B4,C4)
    Kalman4=np.transpose(np.squeeze([B4, A4 @ B4, A4@A4@B4],axis=2))
    print('n=3, rank(k)=')
    print(np.linalg.matrix_rank(Kalman4))
    # -------------------------------- Symulacje ------------------
    t=np.linspace(0,15,1000)
    In1=np.ones_like(t)
    In2=2*np.ones_like(t)
    In3=np.sin(t)-1/2*np.ones_like(t)
    # UKLAD 1 ------- 3 rózne wejscia ------------------
    t1_1,y1_1,x1_1=scipy.signal.lsim(SS1,In1,t)
    t1_2,y1_2,x1_2=scipy.signal.lsim(SS1,In2,t)
    t1_3,y1_3,x1_3=scipy.signal.lsim(SS1,In3,t)
    # UKLAD 2 ------- 3 rózne wejscia ------------------
    t2_1,y2_1,x2_1=scipy.signal.lsim(SS2,In1,t)
    t2_2,y2_2,x2_2=scipy.signal.lsim(SS2,In2,t)
    t2_3,y2_3,x2_3=scipy.signal.lsim(SS2,In3,t)
    # UKLAD 3 ------- 3 rózne wejscia ------------------
    t3_1,y3_1,x3_1=scipy.signal.lsim(SS3,In1,t)
    t3_2,y3_2,x3_2=scipy.signal.lsim(SS3,In2,t)
    t3_3,y3_3,x3_3=scipy.signal.lsim(SS3,In3,t)
    # UKLAD 4 ------- 3 rózne wejscia ------------------
    t4_1,y4_1,x4_1=scipy.signal.lsim(SS4,In1,t)
    t4_2,y4_2,x4_2=scipy.signal.lsim(SS4,In2,t)
    t4_3,y4_3,x4_3=scipy.signal.lsim(SS4,In3,t)
    # -----------------------WYKRESY 1 UKLAD ------------------
    plt.figure(figsize=(8, 9), dpi=130)
    plt.subplot(3,1,1)
    plt.plot(t1_1,x1_1)
    plt.grid()
    plt.legend(['Uc1','Uc2'])
    plt.title('Odpowiedz ukladu na 1H')
    plt.subplot(3,1,2)
    plt.plot(t1_2,x1_2)
    plt.grid()
    plt.legend(['Uc1','Uc2'])
    plt.title('Odpowiedz ukladu na 2*1H')
    plt.subplot(3,1,3)
    plt.plot(t1_3,x1_3)
    plt.grid()
    plt.legend(['Uc1','Uc2'])
    plt.title('Odpowiedz ukladu na sin(t)-1/2')
    plt.show()
    # -----------------------WYKRESY 2 UKLAD ------------------
    plt.figure(figsize=(8, 9), dpi=130)
    plt.subplot(3,1,1)
    plt.plot(t2_1,x2_1)
    plt.grid()
    plt.legend(['Uc1','Uc2','Uc3'])
    plt.title('Odpowiedz ukladu na 1H')
    plt.subplot(3,1,2)
    plt.plot(t2_2,x2_2)
    plt.grid()
    plt.legend(['Uc1','Uc2','Uc3'])
    plt.title('Odpowiedz ukladu na 2*1H')
    plt.subplot(3,1,3)
    plt.plot(t2_3,x2_3)
    plt.grid()
    plt.legend(['Uc1','Uc2','Uc3'])
    plt.title('Odpowiedz ukladu na sin(t)-1/2')
    plt.show()
    # -----------------------WYKRESY 3 UKLAD ------------------
    plt.figure(figsize=(8, 9), dpi=130)
    plt.subplot(3,1,1)
    plt.plot(t3_1,x3_1)
    plt.grid()
    plt.legend(['iL1','Uc1','iL2','Uc2'])
    plt.title('Odpowiedz ukladu na 1H')
    plt.subplot(3,1,2)
    plt.plot(t3_2,x3_2)
    plt.grid()
    plt.legend(['iL1','Uc1','iL2','Uc2'])
    plt.title('Odpowiedz ukladu na 2*1H')
    plt.subplot(3,1,3)
    plt.plot(t3_3,x3_3)
    plt.grid()
    plt.legend(['iL1','Uc1','iL2','Uc2'])
    plt.title('Odpowiedz ukladu na sin(t)-1/2')
    plt.show()
    # -----------------------WYKRESY 4 UKLAD ------------------
    plt.figure(figsize=(8, 9), dpi=130)
    plt.subplot(3,1,1)
    plt.plot(t4_1,x4_1)
    plt.grid()
    plt.legend(['i1','i3','Uc1'])
    plt.title('Odpowiedz ukladu na 1H')
    plt.subplot(3,1,2)
    plt.plot(t4_2,x4_2)
    plt.grid()
    plt.legend(['i1','i3','Uc1'])
    plt.title('Odpowiedz ukladu na 2*1H')
    plt.subplot(3,1,3)
    plt.plot(t4_3,x4_3)
    plt.grid()
    plt.legend(['i1','i3','Uc1'])
    plt.title('Odpowiedz ukladu na sin(t)-1/2')
    plt.show()
    # --------------------- Macierze Kalmara XD -----------------
def Zadanie2():
    # Zadanie 2.1
    # ------------------------------- Uklad 2 ----------------------------------
    C1=1
    C2=2
    C3=3
    #Macierze stanu 2
    A2=np.array([[-1/C1,0,0],[0,-1/C2,0],[0,0,-1/C3]])
    B2=np.array([[1/C1],[1/C2],[1/C3]])
    C2=np.eye(3)
    # Kalman2=[B2,np.dot(A2,B2),np.dot(np.dot(A2,A2),B2)]
    # print("Macierz kalmana powinna mieć rank(K)=3")
    # print(np.linalg.matrix_rank(Kalman2))
    SS2=scipy.signal.StateSpace(A2,B2,C2)
    # ------------------------------- Uklad 4 ----------------------------------
    L1=0.5
    L2=1
    C1=2
    A4=np.array([[-2/L1,0,-1/L1],[0,0,1/L2],[1/C1,-1/C1,-1/C1]])
    B4=np.array([[1/L1],[0],[0]])
    C4=np.eye(3)
    # Kalman4=[B4,np.dot(A4,B4),np.dot(np.dot(A4,A4),B4)]
    # print("Macierz kalmana powinna mieć rank(K)=3")
    # print(np.linalg.matrix_rank(Kalman4))
    SS4=scipy.signal.StateSpace(A4,B4,C4)

    # --------------------------- Wyznaczenie postaci sterowalnej ------------------------
    A2s=np.array([[0,1,0],[0,0,1],[-1/6,-1,-11/6]])
    B2s=np.array([[0],[0],[1]])
    A4s=np.array([[0,1,0],[0,0,1],[-2,-7/2,-9/2]])
    B4s=B2s
    P2_inv=np.hstack([B2, A2 @ B2, A2@A2@B2]) @ np.linalg.inv(np.hstack([B2s,A2s@B2s,A2s@A2s@B2s]))
    P4_inv=np.hstack([B4, A4 @ B4, A4@A4@B4]) @ np.linalg.inv(np.hstack([B2s,A2s@B2s,A2s@A2s@B2s]))
    P2=np.linalg.inv(P2_inv)
    P4=np.linalg.inv(P4_inv)
    A2_ster=P2 @ A2 @ P2_inv
    B2_ster=P2 @ B2
    C2_ster=C2 @ P2_inv
    A4_ster=P4 @ A4 @ P4_inv
    B4_ster=P4 @ B4
    C4_ster=C4 @ P4_inv
    SS2_ster=scipy.signal.StateSpace(A2_ster,B2_ster,C2_ster)
    SS4_ster=scipy.signal.StateSpace(A4_ster,B4_ster,C4_ster)
    # ------------------------- Symulacje czasowe ----------------------------
    t=np.linspace(0,15,1000)
    In1=np.ones_like(t)
    In2=2*np.ones_like(t)
    In3=np.sin(t)*np.ones_like(t)-1/2*np.ones_like(t)
    t2_1,y2_1,x2_1=scipy.signal.lsim(SS2,In1,t)
    t2_2,y2_2,x2_2=scipy.signal.lsim(SS2,In2,t)
    t2_3,y2_3,x2_3=scipy.signal.lsim(SS2,In3,t)
    ts2_1,ys2_1,xs2_1=scipy.signal.lsim(SS2_ster,In1,t)
    ts2_2,ys2_2,xs2_2=scipy.signal.lsim(SS2_ster,In2,t)
    ts2_3,ys2_3,xs2_3=scipy.signal.lsim(SS2_ster,In3,t)

    t4_1,y4_1,x4_1=scipy.signal.lsim(SS4,In1,t)
    t4_2,y4_2,x4_2=scipy.signal.lsim(SS4,In2,t)
    t4_3,y4_3,x4_3=scipy.signal.lsim(SS4,In3,t)
    ts4_1,ys4_1,xs4_1=scipy.signal.lsim(SS4_ster,In1,t)
    ts4_2,ys4_2,xs4_2=scipy.signal.lsim(SS4_ster,In2,t)
    ts4_3,ys4_3,xs4_3=scipy.signal.lsim(SS4_ster,In3,t)
    # -------------------------------- Wykresy dla 2 układu ------------------------
    plt.figure(figsize=(14, 8), dpi=130)
    plt.title('Uklad 2')
    plt.subplot(3,2,1)
    plt.plot(t2_1,x2_1)
    plt.plot(t2_1,y2_1)
    plt.grid()
    plt.legend(['Uc1','Uc2','Uc3','y'])
    plt.title('Odpowiedz ukladu na 1H')
    plt.subplot(3,2,3)
    plt.plot(t2_2,x2_2)
    plt.plot(t2_2,y2_2)
    plt.grid()
    plt.legend(['Uc1','Uc2','Uc3','y'])
    plt.title('Odpowiedz ukladu na 2*1H')
    plt.subplot(3,2,5)
    plt.plot(t2_3,x2_3)
    plt.plot(t2_3,y2_3)
    plt.grid()
    plt.legend(['Uc1','Uc2','Uc3','y'])
    plt.title('Odpowiedz ukladu na sin(t)-1/2')
    plt.subplot(3,2,2)
    plt.plot(ts2_1,xs2_1)
    plt.plot(ts2_1,ys2_1)
    plt.grid()
    plt.legend(['xs1','xs2','xs3','y'])
    plt.title('Odpowiedz ukladu na 1H')
    plt.subplot(3,2,4)
    plt.plot(ts2_2,xs2_2)
    plt.plot(ts2_2,ys2_2)
    plt.grid()
    plt.legend(['xs1','xs2','xs3','y'])
    plt.title('Odpowiedz ukladu na 2*1H')
    plt.subplot(3,2,6)
    plt.plot(ts2_3,xs2_3)
    plt.plot(ts2_3,ys2_3)
    plt.grid()
    plt.legend(['xs1','xs2','xs3','y'])
    plt.title('Odpowiedz ukladu na sin(t)-1/2')
    plt.show()
    #---------------------------------- Wykresy dla 4 układu ---------------------
    plt.figure(figsize=(14, 8), dpi=130)
    plt.title('Uklad 4')
    plt.subplot(3,2,1)
    plt.plot(t4_1,x4_1)
    plt.plot(t4_1,y4_1)
    plt.grid()
    plt.legend(['i1','i3','Uc1','y'])
    plt.title('Odpowiedz ukladu na 1H')
    plt.subplot(3,2,3)
    plt.plot(t4_2,x4_2)
    plt.grid()
    plt.legend(['i1','i3','Uc1','y'])
    plt.title('Odpowiedz ukladu na 2*1H')
    plt.subplot(3,2,5)
    plt.plot(t4_3,x4_3)
    plt.grid()
    plt.legend(['i1','i3','Uc1','y'])
    plt.title('Odpowiedz ukladu na sin(t)-1/2')
    plt.subplot(3,2,2)
    plt.plot(ts4_1,xs4_1)
    plt.grid()
    plt.legend(['xs1','xs2','xs3','y'])
    plt.title('Odpowiedz ukladu na 1H')
    plt.subplot(3,2,4)
    plt.plot(ts4_2,xs4_2)
    plt.grid()
    plt.legend(['xs1','xs2','xs3','y'])
    plt.title('Odpowiedz ukladu na 2*1H')
    plt.subplot(3,2,6)
    plt.plot(ts4_3,xs4_3)
    plt.grid()
    plt.legend(['xs1','xs2','xs3','y'])
    plt.title('Odpowiedz ukladu na sin(t)-1/2')
    plt.show()
    print(SS2_ster)
    print(SS4_ster)

#------------------------------------------ ZADANIE 3 ----------------------------
C1=1
C2=2
C3=3
#Macierze stanu 2
A2=np.array([[-1/C1,0,0],[0,-1/C2,0],[0,0,-1/C3]])
B2=np.array([[1/C1],[1/C2],[1/C3]])
C2=np.array([[1,0,0]])
A2s=np.array([[0,1,0],[0,0,1],[-1/6,-1,-11/6]])
B2s=np.array([[0],[0],[1]])
P2_inv=np.hstack([B2, A2 @ B2, A2@A2@B2]) @ np.linalg.inv(np.hstack([B2s,A2s@B2s,A2s@A2s@B2s]))
P2=np.linalg.inv(P2_inv)
A2_ster=P2 @ A2 @ P2_inv
B2_ster=P2 @ B2
C2_ster=C2 @ P2_inv
SS2_ster=scipy.signal.StateSpace(A2_ster,B2_ster,C2_ster)
K2=scipy.signal.place_poles(A2_ster,B2_ster,[-5,-1,-2])

def state_equations(t, x):
    # Definicja równań stanu
    dxdt = np.dot(A2_ster - np.dot(B2_ster,K2.gain_matrix), x)  # dx/dt = (A-BK)x
    return dxdt

def Zadanie3():

    print(K2.gain_matrix)

    # -------------------- Obliczenia i wykresy ---------------------
    # Symulacja układu
    t=np.linspace(0,15,1000)
    x0=[0.1,0.2,0.3]
    sol = solve_ivp(state_equations, [t[0], t[-1]], x0, t_eval=t)
    
    # Wyniki symulacji
    x_simulated = sol.y
    plt.figure(figsize=(10, 6))
    for i in range(3):
        plt.plot(t, x_simulated[i], label=f'x_{i+1}')

    plt.xlabel('Czas')
    plt.ylabel('Wartość zmiennych stanu')
    plt.legend()
    plt.grid(True)
    plt.title('Przebiegi zmiennych stanu w czasie')
    plt.show()
    

Zadanie2 ()