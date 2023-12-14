close all; clc; 
%Parametry stałe
%Wektor Q_Destined
qd_Vec=[1 2 5];
x0=zeros(5,1); % warunki poczatkowe
Q_vec=[1]; %wektor wartości Q LQR
R_vec=[1]; %wektor wartości R LQR
%Wartości do równan stanu
    R_a = 0.5;
    C = 0.5;
    L = 0.2;
    A = [0 1; -1/(L*C) -R_a/L];
    B = [0; 1/L];
    % Wektory czasu dla róznyych symulacji
    tspan_model=linspace(0,5,1000); %ten mozna tez uzyc do dynamicznego modelu
    tspan_ricatti=linspace(5,0,1000);
for k = 1:length(Q_vec)
    Q=Q_vec(k);
    for j = 1:length(R_vec)
        R=R_vec(j);
        for i =1:3
            q_destined=qd_Vec(i);
            [P, LQR, wart_wlasne] = care(A,B,Q,R);
            global p_ricatti;
            [~,p_ricatti]=ode45(@ricatti,tspan_ricatti,eye(2,2));
            [t,y]=ode45(@(t,y) model(t, y,P, R, Q,q_destined),tspan_model,x0);
            [t_dyn,y_dyn] = ode45(@(t,y) model_dynamic(t,y,R,Q,q_destined),tspan_model,x0);
            %Wartość wskaźnika ---------------------------------------------
            J_opt=y(end,3);
            % Wykres  nieskonczony horyzont -------------------------------------------------------
            figure('Position',[0 0 1920 1080]);
            subplot(3,1,1)
            plot(t,y(:,1:2),'LineWidth',2);
            hold on
            plot(t,y(:,4:5),'LineWidth',2);
            plot(t,y(:,3),'LineWidth',2);

            grid on
            xlabel("czas [s]",'Interpreter','latex')
            ylabel("Wartosci zmiennych stanu",'Interpreter','latex')
            lgd=legend('$x_1$','$x_2$','e','$\dot{e}$','J','','Interpreter','latex');
            title("Model nieskonczony horyzont")
            %Wykres dynamicznego skonczony horyzont
            subplot(3,1,2)
            plot(t_dyn,y_dyn(:,1:2),'LineWidth',2);
            hold on
            plot(t_dyn,y_dyn(:,4:5),'LineWidth',2);
            plot(t_dyn,y_dyn(:,3),'LineWidth',2);
            grid on
            xlabel("czas [s]",'Interpreter','latex')
            ylabel("Wartosci zmiennych stanu",'Interpreter','latex')
            lgd2=legend('$x_1$','$x_2$','e','$\dot{e}$','J','Interpreter','latex');
            title("Model skonczony horyzont")
            sgtitle(['model nieskonczony horyzont Q=' num2str(Q) ' R=' num2str(R) ' q_d=' num2str(q_destined)])
            subplot(3,1,3)
            t_ricatti=linspace(0,5,length(p_ricatti));
            plot(t_ricatti, p_ricatti);
            
        end
    end
end
%%




% Funkcja obliczająca równanie ricattiego w danym zestawie czasowym
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

% Model schematu sterowania
function dxdt = model(t, x, P, R, Q,qd)
    R_a = 0.5;
    C = 0.5;
    L = 0.2;
    A = [0 1; -1/(L*C) -R_a/L];
    B = [0; 1/L];
    K=R^(-1)*B.'*P;
    Axd=[0;-1/(L*C)*qd]; 
    dxdt = zeros(5, 1);
    x_vec=x(1:2);
    e_vec=x(4:5);
    u_e=-K*e_vec;
    u_c = 1/C*qd;
    u = -u_e + u_c ;
    % Zmienne stanu
    dxdt(1:2) = A * x_vec + B * u; 
    % Wskaźnik całkowy 
    dxdt(3) = e_vec.'*Q*e_vec + u_e.'*R*u_e;
    % Wyznaczanie uchybów z nowego rownania
    dxdt(4:5)= A * e_vec  - B * u - Axd;
end

% Funkcja modelu dynamicznego ze zmiennymi LQR
function dxdt = model_dynamic(t, x, R, Q,qd)
    n=1000-cast(200*t,'uint32');
%     disp(['n = ' num2str(n)])
    global p_ricatti;
    if n==0
        P=[p_ricatti(1,1),p_ricatti(1,3);p_ricatti(1,2) p_ricatti(1,4)];
    else
            P=[p_ricatti(n,1),p_ricatti(n,3);p_ricatti(n,2) p_ricatti(n,4)];
    end
    R_a = 0.5;
    C = 0.5;
    L = 0.2;
    A = [0 1; -1/(L*C) -R_a/L];
    B = [0; 1/L];
    K=R^(-1)*B.'*P;
    Axd=[0;-1/(L*C)*qd]; 
    dxdt = zeros(5, 1);
    x_vec=x(1:2);
    e_vec=x(4:5);
    u_e=-K*e_vec;
    u_c = 1/C*qd;
    u = -u_e + u_c ;
    % Zmienne stanu
    dxdt(1:2) = A * x_vec + B * u; 
    % Wskaźnik całkowy 
    dxdt(3) = e_vec.'*Q*e_vec + u_e.'*R*u_e;
    % Wyznaczanie uchybów z nowego rownania
    dxdt(4:5)= A * e_vec  - B * u - Axd;
end