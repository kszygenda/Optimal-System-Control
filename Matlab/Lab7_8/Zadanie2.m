clc; clear all; close all;
% Parametry
    Q=eye(2,2);
    R=0.2;
    R_a=0.5;
    C=0.5;
    L=0.2;
    A = [0 1; -1/(L*C) -R_a/L];
    B = [0 ;1/L];
% Wektory wartości Q i R do testowania wpływu ich wartości 
R_vec=linspace(0.01,1,5);
Q_vec=linspace(0.01,1,5);

    tspan=linspace(0,5,1000);
    x0=[1 1 0];
    fig=figure('Position',[0 0 1920 1080]);
% Podwójna pętla tak aby utworzyć wykres ze wszystkimi odpowiedziami
for j=1:length(Q_vec)
for i=1:length(R_vec)
    %przepisanie wartości do obliczeń
    Q_new=Q_vec(j)*Q;
    R=R_vec(i);
    % Rozwiązanie równania ricattiego
    [P, LQR, wart_wlasne] = care(A,B,Q_new,R);

    % Symulacja model
    [t,y]=ode45(@(t,y) model(t, y,P, R, Q_new),tspan,x0);
    % Wykres
    subplot(5,5,(j-1)*5+i)
    plot(t,y,'LineWidth',2);
    grid on
    xlabel("czas [s]",'Interpreter','latex')
    ylabel("Wartosci zmiennych stanu",'Interpreter','latex')
    lgd=legend('$x_1$','$x_2$','$J_{LQR}$','Interpreter','latex');
    title(['Przebieg dla R=' num2str(R) ' i Q_{new} = ' num2str(Q_vec(j))])
end
end


% Model
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
%% Pytania ze skryptu oraz odpowiedzi
% ---------------------- ZADANIE 2 --------------------------------

% Pytanie:
% Jakie algorytmy są wykorzystywane w funkcji solve_continuous_are?
% Odpowiedź:
% Z dokumentacji można przeczytać, że "The equation is solved by forming
% the extended hamiltonian matrix pencil, as described in [1], given by the
% block matrices" and using QZ decomposition method. 

% Pytanie: 
% Jaki jest charakter odpowiedzi skokowej obiektu?
% Odpowiedź:
% Charakter odpowiedzi skokowej obiektu jest zależny od czasu. Na początku, gdy kondensator jest rozładowany, prąd ma wartości bliskie zeru. Gdy kondensator się ładuje, prąd zaczyna maleć i potem oscylować, aż osiągnie wartość ustaloną.

% Pytanie
% Czy wszystkie zmienne stanu zbiegają do wartości zadanych?
% Odpowiedź
% 

