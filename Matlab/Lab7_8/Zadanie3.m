clc; clear all; close all;
% Parametry
    Q=eye(2,2);
    R=0.2;
    R_a=0.5;
    C=0.5;
    L=0.2;
    A = [0 1; -1/(L*C) -R_a/L];
    B = [0 ;1/L];
% Macierze stanu
t_end_v=[1 2 5];
t_end=0; S=0;
S_vec=[1 100];
for i=1:3
    t_end=t_end_v(i);
    for j=1:2
    S=S_vec(j);
tspan=linspace(t_end,0,1000);
global p_ricatti;
p0=S*eye(2,2);
[t,p_ricatti]=ode45(@ricatti,tspan,p0);
% figure;
% plot(t,p_ricatti)
% legend('p1','p2','p3','p4');
% grid on

%Symulacja modelu dynamicznego.
x0=[1 1 0];
tspan_model=linspace(0,5,1000);
[t,y]=ode45(@(t,y) model_dyn(t, y, R, Q),tspan_model,x0);
siem=figure('Position',[0 0 1920 1080]);
subplot(2,1,1)
plot(t,y(:,1:2),'LineWidth',2);
grid on
hold on
J_opt=y(:,3) + (y(end,1:2)*S*y(end,1:2).');
plot(t,J_opt)
xlabel("czas [s]",'Interpreter','latex')
ylabel("Wartosci zmiennych stanu",'Interpreter','latex')
lgd=legend('$x_1$','$x_2$','$J_{LQR}$','Interpreter','latex');
title('model skonczony horyzont');

%symulacja modelu
[P, LQR, wart_wlasne] = care(A,B,Q,R);
[t,y]=ode45(@(t,y) model(t, y,P, R, Q),tspan_model,x0);
% Wykres
subplot(2,1,2)
plot(t,y,'LineWidth',2);
grid on
xlabel("czas [s]",'Interpreter','latex')
ylabel("Wartosci zmiennych stanu",'Interpreter','latex')
lgd=legend('$x_1$','$x_2$','$J_{LQR}$','Interpreter','latex');
title('model nieskonczony horyzont')
sgtitle(['Wykres dla S = ' num2str(S) 'oraz t_{end} = ' num2str(t_end)])
    end
end






function dxdt = model(t, x, P, R, Q)
    R_a = 0.5;
    C = 0.5;
    L = 0.2;
    A = [0 1; -1/(L*C) -R_a/L];
    B = [0; 1/L];
    K=R^(-1)*B.'*P;
    dxdt = zeros(3, 1);
    x_vec=x(1:2);
    u = -K*x_vec;
    dxdt(1:2) = A * x_vec + B * u;
    dxdt(3) = x_vec.'*Q*x_vec + u.'*R*u; 
end

function dpdt = ricatti(t,P)
    P=reshape(P,[2,2]);
    Q=eye(2,2);
    R=1;
    R_a=0.5;
    C=0.5;
    L=0.2;
    A = [0 1; -1/(L*C) -R_a/L];
    B = [0 ;1/L];
    dpdt=zeros(2,2);
    dpdt =-(P*A - P*B*(1/R)*(B.')*P + (A.'*P) + Q); 
    dpdt=reshape(dpdt,[4,1]);

end

function dxdt = model_dyn(t, x, R, Q)
    n=1000-cast(200*t,'uint32');
%     disp(['n = ' num2str(n)])
    global p_ricatti;
    if n==0
        P=[p_ricatti(1,1),p_ricatti(1,2);p_ricatti(1,3) p_ricatti(1,4)];
    else
            P=[p_ricatti(n,1),p_ricatti(n,2);p_ricatti(n,3) p_ricatti(n,4)];
    end

    R_a = 0.5;
    C = 0.5;
    L = 0.2;
    A = [0 1; -1/(L*C) -R_a/L];
    B = [0; 1/L];
    K=R^(-1)*B.'*P;
    dxdt = zeros(3, 1);
    x_vec=x(1:2);
    u = -K*x_vec;
    dxdt(1:2) = A * x_vec + B * u;
    dxdt(3) = x_vec.'*Q*x_vec + u.'*R*u; 
end
