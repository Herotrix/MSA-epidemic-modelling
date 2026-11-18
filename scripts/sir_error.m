
clear;
clc;

% Locate project folder
script_dir = fileparts(mfilename('fullpath'));
project_dir = fileparts(script_dir);
cd(project_dir);

% Read India second-wave data
data_file = fullfile(project_dir, 'data', 'india_second_wave.csv');
raw_data = csvread(data_file, 1, 1);
observed_cases = raw_data(:, 2);
observed_cases = observed_cases(:);

% Model parameters
N = 1e7;       % Effective population assumption
gamma = 1/7;   % Recovery rate per day
I0 = 1000;
R_initial = 0;
y0 = [N-I0-R_initial; I0; R_initial];

t = (0:length(observed_cases)-1)';
options = odeset('RelTol', 1e-4, 'AbsTol', 1e-2);
weights = 1 ./ max(observed_cases, 1000);

% Search a finite grid of transmission rates
beta_values = 0.05:0.01:0.80;
scale_values = logspace(0, 3, 31);

best_error = Inf;
best_beta = NaN;
best_scale = NaN;

for beta = beta_values
    sir_odes = @(time, y) [
        -beta*y(1)*y(2)/N;
         beta*y(1)*y(2)/N - gamma*y(2);
         gamma*y(2)
    ];

    [~, y] = ode45(sir_odes, t, y0, options);
    incidence = beta .* y(:,1) .* y(:,2) ./ N;

    % Best scale for each beta can be calculated directly
    scale = sum(weights .* incidence .* observed_cases) / ...
            sum(weights .* incidence.^2);
    scale = max(scale, 0);

    predicted = scale .* incidence;
    err = mean(weights .* (predicted - observed_cases).^2);

    if isfinite(err) && err < best_error
        best_error = err;
        best_beta = beta;
        best_scale = scale;
    end
end

if ~isfinite(best_error)
    error('Calibration failed. Check the input CSV and model values.');
end

R0_estimate = best_beta / gamma;

fprintf('Estimated beta: %.5f per day\n', best_beta);
fprintf('Fixed gamma: %.5f per day\n', gamma);
fprintf('Implied beta/gamma: %.3f\n', R0_estimate);
fprintf('Estimated reporting scale: %.5f\n', best_scale);
fprintf('Weighted mean squared error: %.3f\n', best_error);

% Save results
out_file = fullfile(project_dir, 'results', 'sir_fit_parameters.csv');
fid = fopen(out_file, 'w');
fprintf(fid, 'beta,gamma,R0_estimate,reporting_scale,weighted_error\n');
fprintf(fid, '%.8f,%.8f,%.8f,%.8f,%.8f\n', ...
    best_beta, gamma, R0_estimate, best_scale, best_error);
fclose(fid);

fprintf('Parameters saved to: %s\n', out_file);
