clear all; close all;
format longE
tspan = linspace(0, 1, 1000);

% Definicja funkcji celu z dodatkowym parametrem tspan
fun = @(a) integral(@(t) model0000(t, a), 0, 1);

% Początkowe wartości dla parametrów a0, a1, a2, a3
initial_guess = [1.000e+00, 1.000e+00, -4.256e-04, 1.000e+00];

% Ograniczenia dla parametrów
A = [1 0 0 0; 1 1 1 1];
b = [1; 3];

% Opcje optymalizatora
options = optimset('TolFun', 1e-6, 'TolX', 1e-6, 'Algorithm', 'interior-point');
options.MaxIter = 50000;

% Przeprowadzenie optymalizacji
optim = fmincon(fun, initial_guess, [], [], A, b, [], [], [], options);

% Wyświetlenie wyników
disp('Optymalne wartości parametrów:');
disp(optim);

function dxdt = fun_xt(t, a)
    a0 = a(1);
    a1 = a(2);
    a2 = a(3);
    a3 = a(4);
    dxdt = a0 + a1 * t + a2 * t.^2 + a3 * t.^3;
end

function dxdt = fun_xpt(t, a)
    a0 = a(1);
    a1 = a(2);
    a2 = a(3);
    a3 = a(4);
    dxdt = a1 + 2 * a2 .* t + 3 * a3 * t.^2;
end

function J = model(t, a)
    J = 24 * fun_xt(t, a) .* t + 2 * fun_xpt(t, a).^2 - 4 * t;
    disp(J);
end



function [t,y] = problem_dyn(a)
a0=a(1);
a1=a(2);
a2=a(3);
a3=a(4);
tspan=[0 1];
[t,y] = ode45(@(t, y) model(t, y, a0, a1, a2, a3),tspan,0);
end

% Jeżeli wykorzystam linear constraint funkcji fmincon
% dla warunkuw x(0)=1 oraz x(1)=3
% to rownania sa nastepujace
% x(0)=a0 == 1
% x(1)=a0+a1+a2+a3 = 3





% % ------------- CHAT GPT PODEJSCIE 1 ------------------------
% % d) Optymalizacja
% % Ograniczenia dla parametrów
% lb = [-Inf, -Inf, -Inf, -Inf];
% ub = [Inf, Inf, Inf, Inf];
% 
% 
% x0 = [1, 1, 1, 1]; % Przykładowe wartości początkowe
% 
% % Optymalizacja parametrów
% options = optimset('Algorithm', 'interior-point');
% optimal_params = fmincon(@(params) 0, x0, [], [], [], [], lb, ub, @(params) constraints(params), options);
% 
% disp('Znalezione zoptymalizowane parametry:');
% disp(optimal_params);
% 
% % a) Funkcja model
% function output = model(t, a0, a1, a2, a3)
%     x = a0 + a1 * t + a2 * t.^2 + a3 * t.^3;  % x(t)
%     x_prim = a1 + 2 * a2 * t + 3 * a3 * t.^2;  % x'(t)
% 
%     integrand = 24 * x .* t + 2 * x_prim.^2 - 4 * t;
%     output = trapz(t, integrand);
% end
% 
% % b) Funkcja problem_dyn
% function result = problem_dyn(a)
%     % Parametry
%     a0 = a(1);
%     a1 = a(2);
%     a2 = a(3);
%     a3 = a(4);
% 
%     % Horyzont czasowy
%     tspan = [0, 1];
% 
%     % Całkowanie
%     result = model(tspan, a0, a1, a2, a3);
% end

% % c) Funkcja constraints
% function [c, ceq] = constraints(a)
%     % Warunki ograniczające - x(0) = 1 oraz x(1) = 3
%     ceq = [problem_dyn(a) - 3; problem_dyn(a) - 1];
%     c = [];  % Brak dodatkowych ograniczeń nierównościowych
% end











% tspan = [0 1];
% a0_val = 1;
% a1_val = 2;
% a2_val = 3;
% a3_val = 4;
% 
% % Warunki początkowe
% x0 = [1; 3]; % x(0) = 1, x(1) = 3
% 
% % Ustawienie opcji ode45
% options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
% 
% % Rozwiązanie równań różniczkowych
% [t, x] = ode45(@(t, x) ode_equations(t, x, a0_val, a1_val, a2_val, a3_val), tspan, x0, options);
% 
% % Wykres wyniku
% plot(t, x(:, 1), t, x(:, 2));
% xlabel('Czas');
% ylabel('Wartość x');
% legend('x_1', 'x_2');
% title('Rozwiązanie równań różniczkowych');
% 
% 
% function dydt = fun_x(t, a0, a1, a2, a3)
%     dydt = a0 + a1 * t + a2 * t^2 + a3 * t^3;
% end
% 
% function dydt = fun_x_prim(t,a1,a2,a3)
%     dydt =  a1 + 2*a2*t + 3*a3*t^2;
% end
% 
% function dxdt = ode_equations(t, x, a0, a1, a2, a3)
%     fun_x_val = a0 + a1 * t + a2 * t^2 + a3 * t^3;
%     fun_x_prim_val = a1 + 2 * a2 * t + 3 * a3 * t^2;
% 
%     dxdt = 24 * fun_x_val * t + 2 * fun_x_prim_val^2 - 4 * t;
% end












