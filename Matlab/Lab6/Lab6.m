clc; clear all; close all;

tspan = [0 5];
x0 = [0 0 0 0];
[t, y] = ode45(@model, tspan, x0);
yd=3;
Zadanie2 = figure('Position', [10 10 1200 800]);
plot(t, y(:, 2), 'LineWidth', 2)
hold on
plot(t, yd-y(:,2), 'LineWidth', 2)
xlabel('czas [s]', 'FontSize', 12)
ylabel('Prąd [A]', 'FontSize', 12)
lgd = legend('y(t)', 'e(t)', 'Interpreter', 'latex');
title('Zadanie 2.1', 'FontSize', 19)
fontsize(lgd, 14, 'points')
grid on

%Zadanie 3.1
T_osc=2.6;
ku=25;
kp_P=0.5*ku; % 12.5
kp_PI=0.45*ku; % 11.2500
ki_PI=0.54*ku*T_osc^(-1); %5.1923
kp_PD=0.8*ku; % 20
kd_PD=0.1*ku*T_osc; % 6.5000
kp_PID=0.6*ku; % 15
ki_PID=1.2*ku*T_osc^(-1); % 11.5385
kd_PID=0.075*ku*T_osc; % 4.8750
function dxdt = model(t, x)
    % Parametry układu
    R1 = 0.2;
    R2 = 5;
    C1 = 0.5;
    L1 = 2;
    L2 = 0.5;
    
    % Regulator PID
    kp = 10;
    ki = kp/0.5;
    kd = 0;
    yd = 3;
    e=yd-x(2);
    e_prim=0 * x(1) - R2/L2 * x(2) + 1/L2 * x(3);
    u=kp*e+ki*x(4)+kd*(-e_prim);
    % Obliczenia
    dxdt = zeros(8, 1);
    
    dxdt(1) = -R1/L1 * x(1) + 0 * x(2) - 1/L1 * x(3) + 1/L1 * u;
    dxdt(2) = 0 * x(1) - R2/L2 * x(2) + 1/L2 * x(3);
    dxdt(3) = 1/C1 * x(1) - 1/C1 * x(2) - 0 * x(3);
    dxdt(4) = e;
    %Kryteria całkowe, po kolei
    dxdt(5) = e^2; %ISE
    dxdt(6) = t*e^2; %ITSE
    dxdt(7) = abs(e); %IAE
    dxdt(8) = t*abs(e); %ITAE
end
%Pytanie 1
%Jak zwiększanie/zmniejszanie wartości każdego z parametrów regulatora PID wpływają na odpowiedź skokową układu?
% No standardowo pidzik nie
%Pytanie 2
%Czy możliwe jest uzyskanie zerowego uchybu ustalonego dla regulatora typu P?
% Nie, granica stabilnosci mowi kaboom
% Pytanie 3
% kp=1 ki=
%Zadanie 3.1


