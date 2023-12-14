import numpy as np
import matplotlib . pyplot as plt
from scipy . integrate import solve_ivp
import scipy.signal
import math

def model1(t,y):
    return t**2
def Zadanie1():
    t=np.linspace(0,10,100)
    y0=0
    # Solving the ODE // model,tspan, initial conditions, time points. Used method is Runge-Kutta 45
    sol=solve_ivp(model1,[0,10],[0],t_eval=t)
    new_y=1/3*t**3 # Analytical solution


    # Plotting the solution
    plt.figure(figsize=(8, 9), dpi=130)
    plt.plot(sol.t,sol.y[0],'b-',label='Numerical solution') # Updated line
    plt.plot(t,new_y,'r--',label='Analytical solution')
    plt.xlabel('Time')
    plt.ylabel('y(t)')
    plt.legend()
    plt.grid()
    plt.title('Solution of ODE')
    plt.show()

def model2(t,x):
    kp=4
    w=4
    psi=0.25
    u=1
    x1 = x[0] #y
    x2 = x[1] #dydt
    dydt2 = (-2*psi/w)*x2 - (1/w)*np.sqrt(x1) + (kp/(w**2))*u 
    return np.array([x2, dydt2])

def Zadanie2():
    t=np.linspace(0,50,1000)
    y0=[0, 0]  # Update initial conditions to include both x1 and x2

    # Solving the ODE // model,tspan, initial conditions, time points. Used method is Runge-Kutta 45
    sol=solve_ivp(model2,[0,50],y0,t_eval=t)  # Update initial conditions argument
    
    # Plotting the solution
    plt.figure(figsize=(8, 9), dpi=130)
    plt.plot(sol.t,sol.y[0],'b-',label='x1 = y') # Updated line
    plt.plot(sol.t,sol.y[1],'r-',label='x2 = dydt') # Updated line
    plt.xlabel('Time [s]')
    plt.ylabel('wartosci x1 i x2')
    plt.legend()
    plt.grid()
    plt.title('Solution of ODE')
    plt.show()

def feedback(t,y,xd):
    # Rozpiska rownania
    e=xd-y[0]
    kp=2
    kob=4
    T=2
    u=kp*e
    uc=np.clip(u,-0.1,0.1)
    uc=u
    dydt=(kob/T)*uc-(1/T)*y[0] # CO DO CHUJA
    return np.array([dydt])

def Zadanie3():
    t=np.linspace(0,20,1000)
    y0=[0]  # Update initial conditions to include both x1 and x2
    xd_vec=[1,2,3]
    # Solving the ODE // model,tspan, initial conditions, time points. Used method is Runge-Kutta 45
    sol1=solve_ivp(feedback,[0,20],y0,t_eval=t,args=(xd_vec[0],))  
    sol2=solve_ivp(feedback,[0,20],y0,t_eval=t,args=(xd_vec[1],)) 
    sol3=solve_ivp(feedback,[0,20],y0,t_eval=t,args=(xd_vec[2],)) 
    # Plotting the solution
    plt.figure(figsize=(8, 9), dpi=130)
    plt.subplot(3,1,1)
    plt.plot(sol1.t,sol1.y[0],'b-',label='x1 = y') # Updated line
    plt.title('xd = 1')
    plt.ylabel('y(t)')
    plt.legend()
    plt.grid()
    plt.subplot(3,1,2)
    plt.plot(sol2.t,sol2.y[0],'r-',label='x1 = y') # Updated line
    plt.title('xd = 2')
    plt.ylabel('y(t)')
    plt.legend()
    plt.grid()
    plt.subplot(3,1,3)
    plt.plot(sol3.t,sol3.y[0],'g-',label='x1 = y') # Updated line
    plt.title('xd = 3')
    plt.ylabel('y(t)')
    plt.legend()
    plt.grid()
    plt.xlabel('Time [s]')
    plt.show()


Zadanie3()