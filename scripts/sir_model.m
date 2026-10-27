
% SIR Epidemic Model
% First simulation for the Modelling and Simulations project

clear;
clc;
close all;

% Population and initial conditions
N = 1000000;      % Total population
I0 = 100;         % Initially infected people
R0 = 0;           % Initially recovered people
S0 = N - I0 - R0; % Initially susceptible people

% Model parameters (illustrative values)
beta = 0.30;      % Transmission rate per day
gamma = 0.10;     % Recovery rate per day

% Simulation time in days
tspan = [0 160];

% Initial state vector: [Susceptible, Infected, Recovered]
y0 = [S0; I0; R0];

% SIR differential equations
sir_odes = @(t, y) [
    -beta * y(1) * y(2) / N;
     beta * y(1) * y(2) / N - gamma * y(2);
     gamma * y(2)
];

% Solve the differential equations
[t, y] = ode45(sir_odes, tspan, y0);

% Extract each population
S = y(:, 1);
I = y(:, 2);
R = y(:, 3);

% Plot the simulation
figure;
plot(t, S, 'b-', 'LineWidth', 2);
hold on;
plot(t, I, 'r-', 'LineWidth', 2);
plot(t, R, 'g-', 'LineWidth', 2);
hold off;

xlabel('Time (days)');
ylabel('Number of people');
title('SIR Epidemic Simulation');
legend('Susceptible', 'Infected', 'Recovered', 'Location', 'best');
grid on;

% Check conservation of population
population_error = max(abs(S + I + R - N));

fprintf('Simulation completed successfully.\n');
fprintf('Peak infected population: %.0f\n', max(I));
fprintf('Day of peak infection: %.1f\n', t(find(I == max(I), 1)));
fprintf('Maximum population conservation error: %.6f\n', population_error);

% Save numerical results
results = [t, S, I, R];

output_file = fullfile(fileparts(pwd()), 'results', 'sir_simulation.csv');
fid = fopen(output_file, 'w');

if fid == -1
    error('Could not open output file: %s', output_file);
end

fprintf(fid, 'day,susceptible,infected,recovered\n');
fprintf(fid, '%.6f,%.6f,%.6f,%.6f\n', results');
fclose(fid);

fprintf('Results saved to: %s\n', output_file);
