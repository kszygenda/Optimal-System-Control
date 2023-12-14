% Define the ODE function
odefun = @(t, y) -y;

% Set the time span
tspan = [0 5];

% Set the initial conditions
y0 = 1;

% Define the event condition (stop when y crosses zero)
eventfun = @(t, y) y;

% Set the options with the event function
options = odeset('Events', eventfun);

% Call ode45 with events
[t, y, te, ye] = ode45(odefun, tspan, y0, options);

% Plot the solution and event points
plot(t, y, 'LineWidth', 2)
hold on
plot(te, ye, 'ro')  % Mark event points in red
xlabel('Time')
ylabel('Solution')
legend('Solution', 'Event Points')
grid on
