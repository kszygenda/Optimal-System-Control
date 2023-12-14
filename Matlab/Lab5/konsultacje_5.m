if Zadanie==3
main()
end
% --------------------- FUNKCJE DO SKRYPTU ------------------------------
function main()
    % Warunki początkowe
    a0_0 = 0;
    a1_0 = 0;
    a2_0 = 0;
    a3_0 = 0;
    
    % dziedziny parametrów
    lb = [-Inf, -Inf, -Inf, -Inf]; 
    ub = [Inf, Inf, Inf, Inf];    

    % Warunki początkowe
    options = optimset('fmincon');
    options.Display = 'iter'; % Wyświetlanie wyników optymalizacji
    a_opt = fmincon(@(a) cost_function(a), [a0_0, a1_0, a2_0, a3_0], [], [], [], [], lb, ub, [], options);
    % wyniki
    disp('Optymalne wartości parametrów:');
    disp(a_opt);
    % Całkowanie wyjścia
    tspan = [0, 1]; % Zakres czasu do ode45
    y0 = 0;         % Warunki początkowe funkcji
    [t_output, y_output] = ode45(@(t, y) problem_dyn(t, y, a_opt), tspan, y0);
    
    % Wykres wyników całkowania
    figure;
    plot(t_output, y_output);
    xlabel('t');
    ylabel('output');
    title('model');
end

function dydt = model(t, a0, a1, a2, a3)
    dydt = a0 + a1 * t + a2 * t.^2 + a3 * t.^3;
end

function cost = cost_function(a)
    % Warunki początkowe - dydt(0) = 1 i dydt(1) = 3
    dydt_0 = model(0, a(1), a(2), a(3), a(4));
    dydt_1 = model(1, a(1), a(2), a(3), a(4));
    
    % Funkcja celu - kwadrat różnicy między warunkami początkowymi a 1 i 3
    cost = (dydt_0 - 1)^2 + (dydt_1 - 3)^2;
end

function dydt = problem_dyn(t, y, a)
    dydt = model(t, a(1), a(2), a(3), a(4));
end