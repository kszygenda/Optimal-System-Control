import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from sympy import *

# czas
t=np.linspace(0, 10, 1000)

# zmienne symboliczne
l = symbols('λ')
k1 = symbols('k1')
k2 = symbols('k2')

# Zadanie 1 ###########################################
# zmienne
C=0.01
Vin=8
R=5
L=0.1

# macierze do zadania 1
A=np.array([[0, 1], [-1/(C*L), -1/(C*R)]])
B=np.array([[0], [Vin/(C*L)]])
C=np.array([1, 0])
# wymusznie do zadania 1
mi=0.5

def model1(x, t):
    return [x[1], A[1, 0]*x[0]+A[1, 1]*x[1]+B[1, 0]*mi, 0, 0]

def zadanie1():
    print('Macierze stanu')
    print(A)
    print(B)
    print(C)
    print('')

    print('Wartości własne macierzy A')
    print(np.linalg.eig(A))
    print('')
    print('Układ stabilny(albo blisko granicy stabilności)\n'
          'Układ będzie miał oscylacje, urojone wartości własne(czyli bieguny transmitancji)\n')
# Zadanie 1 ###########################################

# Zadanie 2 ###########################################
def zadanie2():
    print('Wyznaczenie macierzy obserwowalności\n')
    R=np.array([[B[0, 0], A[0, 0]*B[0, 0]+A[0, 1]*B[1, 0]], [B[1, 0], A[1, 0]*B[0, 0]+A[1, 1]*B[1, 0]]])
    print('Macierz sterowalności układu\n')
    print(R)
    print('Rząd macierzy\n')
    print(np.linalg.matrix_rank(R))
    print('Rząd macierzy sterowalności równy rzędowi układu, układ sterowalny\n')
# Zadanie 2 ###########################################

# Zadanie 3 ###########################################
def wyznaczWzmocnienia(param):
    print('Równanie charakterystyczne')
    print('λ**2', '+', 'λ*(8000K2+20)', '+', '8000K1+1000')
    print()
    if param == 1:
        print('(λ+omega0) do kwadratu <------------')
        print(l ** 2 + 2 * l + 1)
        print('Rozwiązania:', 'K1=-0.124875', 'K2=-0.00225')
    if param == 5:
        print('(λ+omega0) do kwadratu <------------')
        print(l ** 2 + 10 * l + 25)
        print('Rozwiązania:', 'K1=-0.121875', 'K2=-0.00125')
    if param == 10:
        print('(λ+omega0) do kwadratu <------------')
        print(l ** 2 + 20 * l + 100)
        print('Rozwiązania:', 'K1=-0.1125', 'K2=0')
    if param == -1:
        print('(λ+omega0) do kwadratu <------------')
        print(l ** 2 - 2 * l + 1)
        print('Rozwiązania:', 'K1=-0.124875', 'K2=-0.00275')
    if param == 2:
        print('(λ+omega0) do kwadratu <------------')
        print(l ** 2 + 4 * l + 4)
        print('Rozwiązania:', 'K1=-0.1245', 'K2=-0.002')

    print('')
    print('Wyznaczenie uD')
    print('Lewa strona równania')
    print(B*symbols('uD'))
    print('Prawa strona równania')
    print(-A@np.array([[5], [0]]))
    print('Po przyrównaniu uD=0.625')

# 4 funkcje model do pokazania tamtych tych i tamtych tam
def model3_1(x, t):
    # Wektor xd = [5, 0]
    # Wzmocnienia do zmieniania
    K1=-0.124875
    K2=-0.00225
    # Przyjęte wektory xd
    e1=5-x[0]
    e2=0-x[1]
    # Wymuszenia
    uD=0.625
    uFB=K1*e1+K2*e2
    u=uFB+uD

    Haha = np.array([[A[0, 0] - B[0, 0] * K1, A[0, 1] - B[0, 0] * K2], [A[1, 0] - B[1, 0] * K1, A[1, 1] - B[1, 0] * K2]])

    return [x[1], A[1, 0]*x[0]+A[1, 1]*x[1]+B[1, 0]*u,
            Haha[0, 0]*e1+Haha[0, 1]*e2, Haha[1, 0]*e1+Haha[1, 1]*e2]

def model3_2(x, t):
    # Wektor xd = [5, 0]
    # Wzmocnienia do zmieniania
    K1=-0.121875
    K2=-0.00125
    # Przyjęte wektory xd
    e1=5-x[0]
    e2=0-x[1]
    # Wymuszenia
    uD=0.625
    uFB=K1*e1+K2*e2
    u=uFB+uD

    Haha = np.array([[A[0, 0] - B[0, 0] * K1, A[0, 1] - B[0, 0] * K2], [A[1, 0] - B[1, 0] * K1, A[1, 1] - B[1, 0] * K2]])

    return [x[1], A[1, 0]*x[0]+A[1, 1]*x[1]+B[1, 0]*u,
            Haha[0, 0]*e1+Haha[0, 1]*e2, Haha[1, 0]*e1+Haha[1, 1]*e2]

def model3_3(x, t):
    # Wektor xd = [5, 0]
    # Wzmocnienia do zmieniania
    K1=-0.1125
    K2=0
    # Przyjęte wektory xd
    e1=5-x[0]
    e2=0-x[1]
    # Wymuszenia
    uD=0.625
    uFB=K1*e1+K2*e2
    u=uFB+uD

    Haha = np.array([[A[0, 0] - B[0, 0] * K1, A[0, 1] - B[0, 0] * K2], [A[1, 0] - B[1, 0] * K1, A[1, 1] - B[1, 0] * K2]])

    return [x[1], A[1, 0]*x[0]+A[1, 1]*x[1]+B[1, 0]*u,
            Haha[0, 0]*e1+Haha[0, 1]*e2, Haha[1, 0]*e1+Haha[1, 1]*e2]

def model3_4(x, t):
    # Wektor xd = [5, 0]
    # Wzmocnienia do zmieniania
    K1=-0.124875
    K2=-0.00275
    # Przyjęte wektory xd
    e1=5-x[0]
    e2=0-x[1]
    # Wymuszenia
    uD=0.625
    uFB=K1*e1+K2*e2
    u=uFB+uD

    Haha = np.array([[A[0, 0] - B[0, 0] * K1, A[0, 1] - B[0, 0] * K2], [A[1, 0] - B[1, 0] * K1, A[1, 1] - B[1, 0] * K2]])

    return [x[1], A[1, 0]*x[0]+A[1, 1]*x[1]+B[1, 0]*u,
            Haha[0, 0]*e1+Haha[0, 1]*e2, Haha[1, 0]*e1+Haha[1, 1]*e2]


# 4 funkcje model do pokazania tych tych tych

def zadanie3():
    print('Wyznaczenie macierzy H')
    H=np.array([[A[0][0] - B[0][0]*k1, A[0][1] - B[0][0]*k2],[A[1][0] - B[1][0]*k1, A[1][1] - B[1][0]*k2]])
    print(H)
    print('')
    # DO WYZNACZENIA WZMOCNIEŃ
    wyznaczWzmocnienia(5)

    sol1=odeint(model3_1, [2, 0, 3, 0], t)

    plt.figure()
    plt.plot(t, sol1[:, 0], 'r')
    plt.plot(t, sol1[:, 1], 'b')
    plt.plot(t, sol1[:, 2], 'c')
    plt.plot(t, sol1[:, 3], 'm')
    plt.xlabel('Czas')
    plt.ylabel('Wartości zmiennych stanu')
    plt.title('Przebiegi dla ω=1')
    plt.legend(labels=['x1', 'x2', 'e1', 'e2'])
    plt.grid()

    sol2=odeint(model3_2, [2, 0, 3, 0], t)

    plt.figure()
    plt.plot(t, sol2[:, 0], 'r')
    plt.plot(t, sol2[:, 1], 'b')
    plt.plot(t, sol2[:, 2], 'c')
    plt.plot(t, sol2[:, 3], 'm')
    plt.xlabel('Czas')
    plt.ylabel('Wartości zmiennych stanu')
    plt.title('Przebiegi dla ω=5')
    plt.legend(labels=['x1', 'x2', 'e1', 'e2'])
    plt.grid()

    sol3=odeint(model3_3, [2, 0, 3, 0], t)

    plt.figure()
    plt.plot(t, sol3[:, 0], 'r')
    plt.plot(t, sol3[:, 1], 'b')
    plt.plot(t, sol3[:, 2], 'c')
    plt.plot(t, sol3[:, 3], 'm')
    plt.xlabel('Czas')
    plt.ylabel('Wartości zmiennych stanu')
    plt.title('Przebiegi dla ω=10')
    plt.legend(labels=['x1', 'x2', 'e1', 'e2'])
    plt.grid()

    sol4=odeint(model3_4, [2, 0, 3, 0], t)

    plt.figure()
    plt.plot(t, sol4[:, 0], 'r')
    plt.plot(t, sol4[:, 1], 'b')
    plt.plot(t, sol4[:, 2], 'c')
    plt.plot(t, sol4[:, 3], 'm')
    plt.xlabel('Czas')
    plt.ylabel('Wartości zmiennych stanu')
    plt.title('Przebiegi dla ω=-1')
    plt.legend(labels=['x1', 'x2', 'e1', 'e2'])
    plt.grid()

    plt.show()
# Zadanie 3 ###########################################

# Zadanie 4 ###########################################
# Szum na wymuszenie u

def model4_1(x, t):
    # Szum na wymuszenie u
    u_szum = 0.5 * np.sin(0.2 * 3.14 * t) ** 2
    # Szum na uchyby e
    e_szum = np.sin(4*3.14*t)**2
    ##########################################
    # Wektor xd = [5, 0]
    # Wzmocnienia do zmieniania
    K1=-0.124875
    K2=-0.00225
    # Przyjęte wektory xd
    e1=5-x[0]
    e2=0-x[1]
    # Wymuszenia
    uD=0.625
    uFB=K1*e1+K2*e2
    u=uFB+uD

    Haha = np.array([[A[0, 0] - B[0, 0] * K1, A[0, 1] - B[0, 0] * K2], [A[1, 0] - B[1, 0] * K1, A[1, 1] - B[1, 0] * K2]])

    # Z szumem na u
    return [x[1], A[1, 0]*x[0]+A[1, 1]*x[1]+B[1, 0]*(u+u_szum)]

    # Z szumem na e
    # Wymuszenia
    # uD=0.625
    # uFB=K1*(e1+e_szum)+K2*(e2+e_szum)
    # u=uFB+uD
    # return [x[1], A[1, 0]*x[0]+A[1, 1]*x[1]+B[1, 0]*u]

def model4_2(x, t):
    # Szum na wymuszenie u
    u_szum = 0.5 * np.sin(0.2 * 3.14 * t) ** 2
    # Szum na uchyby e
    e_szum = np.sin(4*3.14*t)**2
    ##########################################
    # Wektor xd = [5, 0]
    # Wzmocnienia do zmieniania
    K1=-0.1245
    K2=-0.002
    # Przyjęte wektory xd
    e1=0-x[0]
    e2=0-x[1]
    # Wymuszenia
    uD=0.625
    uFB=K1*e1+K2*e2
    u=uFB+uD

    Haha = np.array([[A[0, 0] - B[0, 0] * K1, A[0, 1] - B[0, 0] * K2], [A[1, 0] - B[1, 0] * K1, A[1, 1] - B[1, 0] * K2]])

    # Z szumem na u
    return [x[1], A[1, 0]*x[0]+A[1, 1]*x[1]+B[1, 0]*(u+u_szum)]
    # Z szumem na e
    # # Wymuszenia
    # uD=0.625
    # uFB=K1*(e1+e_szum)+K2*(e2+e_szum)
    # u=uFB+uD
    # return [x[1], A[1, 0]*x[0]+A[1, 1]*x[1]+B[1, 0]*u]

def model4_3(x, t):
    # Szum na wymuszenie u
    u_szum = 0.5 * np.sin(0.2 * 3.14 * t) ** 2
    # Szum na uchyby e
    e_szum = np.sin(4*3.14*t)**2
    ##########################################
    # Wektor xd = [5, 0]
    # Wzmocnienia do zmieniania
    K1=-0.121875
    K2=-0.00125
    # Przyjęte wektory xd
    e1=0-x[0]
    e2=0-x[1]
    # Wymuszenia
    uD=0.625
    uFB=K1*e1+K2*e2
    u=uFB+uD

    Haha = np.array([[A[0, 0] - B[0, 0] * K1, A[0, 1] - B[0, 0] * K2], [A[1, 0] - B[1, 0] * K1, A[1, 1] - B[1, 0] * K2]])

    # Z szumem na u
    return [x[1], A[1, 0]*x[0]+A[1, 1]*x[1]+B[1, 0]*(u+u_szum)]

    # Z szumem na e
    # # Wymuszenia
    # uD=0.625
    # uFB=K1*(e1+e_szum)+K2*(e2+e_szum)
    # u=uFB+uD
    # return [x[1], A[1, 0]*x[0]+A[1, 1]*x[1]+B[1, 0]*u]

def zadanie4():

    sol1=odeint(model4_1, [0, 0], t)

    plt.figure()
    plt.plot(t, sol1[:, 0], 'r')
    plt.plot(t, sol1[:, 1], 'b')
    plt.plot(t, 5-sol1[:, 0], 'c')
    plt.plot(t, 0-sol1[:, 1], 'm')
    plt.xlabel('Czas')
    plt.ylabel('Wartości zmiennych stanu')
    plt.title('Przebiegi dla ω=1')
    plt.legend(labels=['x1', 'x2', 'e1', 'e2'])
    plt.grid()

    sol2=odeint(model4_2, [0, 0], t)

    plt.figure()
    plt.plot(t, sol2[:, 0], 'r')
    plt.plot(t, sol2[:, 1], 'b')
    plt.plot(t, 5-sol2[:, 0], 'c')
    plt.plot(t, 0-sol2[:, 1], 'm')
    plt.xlabel('Czas')
    plt.ylabel('Wartości zmiennych stanu')
    plt.title('Przebiegi dla ω=2')
    plt.legend(labels=['x1', 'x2', 'e1', 'e2'])
    plt.grid()

    sol3=odeint(model4_3, [0, 0], t)

    plt.figure()
    plt.plot(t, sol3[:, 0], 'r')
    plt.plot(t, sol3[:, 1], 'b')
    plt.plot(t, 5-sol3[:, 0], 'c')
    plt.plot(t, 0-sol3[:, 1], 'm')
    plt.xlabel('Czas')
    plt.ylabel('Wartości zmiennych stanu')
    plt.title('Przebiegi dla ω=5')
    plt.legend(labels=['x1', 'x2', 'e1', 'e2'])
    plt.grid()

    plt.show()
# Zadanie 4 ###########################################

if __name__ == '__main__':
    zadanie1()
    zadanie2()
    zadanie3()
    zadanie4()
