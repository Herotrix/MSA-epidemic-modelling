
clear;
clc;
close all;

% Locate project folder
script_dir = fileparts(mfilename('fullpath'));
project_dir = fileparts(script_dir);
cd(project_dir);

% Load observed daily cases
data_file = fullfile(project_dir, 'data', 'india_second_wave.csv');
fid = fopen(data_file, 'r');

if fid == -1
    error('Could not open the India second-wave dataset.');
end

header = fgetl(fid);
raw_data = textscan(fid, '%q%q%q', 'Delimiter', ',');
fclose(fid);

observed_cases = str2double(raw_data{3});
observed_cases = observed_cases(:);

n = length(observed_cases);
t = (0:n-1)';

% Model assumptions
N = 1e7;
gamma = 1/7;
I0 = 1000;
R_initial = 0;

y0 = [N-I0-R_initial; I0; R_initial];

options = odeset('RelTol', 1e-4, 'AbsTol', 1e-2);
weights = 1 ./ max(observed_cases, 1000);

% Search ranges
beta1_values = 0.12:0.04:0.48;
beta2_values = 0.04:0.04:0.32;
change_days = 15:5:100;

best_error = Inf;
best_beta1 = NaN;
best_beta2 = NaN;
best_change_day = NaN;
best_scale = NaN;
best_prediction = [];

% Grid search
for beta1 = beta1_values
    for beta2 = beta2_values
        for change_day = change_days

            k = change_day + 1;

            % Before transmission-rate change
            rhs1 = @(time, y) [
                -beta1*y(1)*y(2)/N;
                 beta1*y(1)*y(2)/N - gamma*y(2);
                 gamma*y(2)
            ];

            [~, y1] = ode45(rhs1, t(1:k), y0, options);

            % After transmission-rate change
            rhs2 = @(time, y) [
                -beta2*y(1)*y(2)/N;
                 beta2*y(1)*y(2)/N - gamma*y(2);
                 gamma*y(2)
            ];

            [~, y2] = ode45(rhs2, t(k:end), y1(end,:)', options);

            % Join both simulation periods
            y = [y1; y2(2:end,:)];

            % Calculate daily new infections
            incidence = zeros(n, 1);

            incidence(1:k) = ...
                beta1 .* y(1:k,1) .* y(1:k,2) ./ N;

            incidence(k+1:end) = ...
                beta2 .* y(k+1:end,1) .* y(k+1:end,2) ./ N;

            % Find best reporting scale
            scale = sum(weights .* incidence .* observed_cases) / ...
                    sum(weights .* incidence.^2);

            scale = max(scale, 0);

            predicted = scale .* incidence;

            % Weighted error
            err = mean(weights .* (predicted - observed_cases).^2);

            if isfinite(err) && err < best_error
                best_error = err;
                best_beta1 = beta1;
                best_beta2 = beta2;
                best_change_day = change_day;
                best_scale = scale;
                best_prediction = predicted;
            end
        end
    end
end

% Display results
fprintf('\nVARIABLE-BETA SIR RESULTS\n');
fprintf('-------------------------\n');
fprintf('Beta before change: %.4f\n', best_beta1);
fprintf('Beta after change:  %.4f\n', best_beta2);
fprintf('Recovery rate:      %.4f\n', gamma);
fprintf('Change day:         %d\n', best_change_day);
fprintf('R0 before change:   %.4f\n', best_beta1/gamma);
fprintf('R0 after change:    %.4f\n', best_beta2/gamma);
fprintf('Reporting scale:    %.4f\n', best_scale);
fprintf('Weighted error:     %.4f\n', best_error);

% Save fitted parameters
out_file = fullfile(project_dir, 'results', ...
                    'sir_variable_beta_metrics.csv');

fid = fopen(out_file, 'w');
fprintf(fid, 'beta_before,beta_after,gamma,change_day,R0_before,R0_after,reporting_scale,weighted_error\n');

fprintf(fid, '%.8f,%.8f,%.8f,%d,%.8f,%.8f,%.8f,%.8f\n', ...
    best_beta1, best_beta2, gamma, best_change_day, ...
    best_beta1/gamma, best_beta2/gamma, best_scale, best_error);

fclose(fid);

% Plot observed and predicted cases
figure;

plot(t, observed_cases, 'k-', 'LineWidth', 2);
hold on;

plot(t, best_prediction, 'r--', 'LineWidth', 2);

xline(best_change_day, 'b:', 'LineWidth', 1.5);

xlabel('Days from 1 March 2021');
ylabel('Daily new cases');
title('India Second Wave: Variable-Beta SIR Fit');

legend('Observed cases', 'Variable-beta SIR', ...
       'Transmission change', 'Location', 'best');

grid on;

% Save plot
print(fullfile(project_dir, 'results', ...
      'sir_variable_beta.png'), '-dpng');

fprintf('\nResults and plot saved successfully.\n');
