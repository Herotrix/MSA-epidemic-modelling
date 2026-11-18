
clear;
clc;
close all;

% Locate project folder
script_dir = fileparts(mfilename('fullpath'));
project_dir = fileparts(script_dir);
cd(project_dir);

% Load calibrated variable-beta parameters
param_file = fullfile(project_dir, 'results', ...
                      'sir_variable_beta_metrics.csv');

params = csvread(param_file, 1, 0);

beta_before = params(1, 1);
beta_after = params(1, 2);
gamma = params(1, 3);
change_day = params(1, 4);
reporting_scale = params(1, 7);

% Model assumptions (same as calibration)
N = 1e7;
I0 = 1000;
R_initial = 0;
S0 = N - I0 - R_initial;
y0 = [S0; I0; R_initial];

options = odeset('RelTol', 1e-4, 'AbsTol', 1e-2);

% Simulation settings
simulation_days = 160;

intervention_strengths = [0, 0.20, 0.40, 0.60];
intervention_start_days = [15, 30, 45, 60, 75];

number_of_scenarios = ...
    length(intervention_strengths) * ...
    length(intervention_start_days);

% Results columns:
% strength, start_day, peak_cases, peak_day,
% cumulative_cases, peak_reduction_percent,
% cumulative_reduction_percent

results = zeros(number_of_scenarios, 7);

scenario = 0;

% Store baseline results for comparison
baseline_peak = NaN;
baseline_cumulative = NaN;

% Run all intervention scenarios
for s = 1:length(intervention_strengths)

    strength = intervention_strengths(s);

    for j = 1:length(intervention_start_days)

        start_day = intervention_start_days(j);

        scenario = scenario + 1;

        % Initial conditions
        y = zeros(simulation_days + 1, 3);
        y(1,:) = y0';

        daily_incidence = zeros(simulation_days, 1);

        % Simulate one day at a time.
        % Beta remains constant within each daily interval.
        for d = 0:simulation_days-1

            % Original fitted transmission rate
            if d < change_day
                beta = beta_before;
            else
                beta = beta_after;
            end

            % Apply additional intervention
            if d >= start_day
                beta = beta * (1 - strength);
            end

            % SIR differential equations
            sir_odes = @(time, state) [
                -beta * state(1) * state(2) / N;
                 beta * state(1) * state(2) / N ...
                    - gamma * state(2);
                 gamma * state(2)
            ];

            % Integrate this one-day interval
            [~, segment] = ode45( ...
                sir_odes, [d, d+1], y(d+1,:)', options);

            y(d+2,:) = segment(end,:);

            % New infections during this daily interval
            daily_incidence(d+1) = ...
                beta * y(d+1,1) * y(d+1,2) / N;
        end

        % Convert model incidence to reported cases
        reported_cases = reporting_scale * daily_incidence;

        % Scenario measurements
        [peak_cases, peak_index] = max(reported_cases);
        peak_day = peak_index - 1;

        cumulative_cases = sum(reported_cases);

        % First scenario is the 0% intervention baseline
        if scenario == 1
            baseline_peak = peak_cases;
            baseline_cumulative = cumulative_cases;
        end

        % Percentage reductions relative to baseline
        peak_reduction = ...
            100 * (baseline_peak - peak_cases) / baseline_peak;

        cumulative_reduction = ...
            100 * (baseline_cumulative - cumulative_cases) ...
            / baseline_cumulative;

        % Store results
        results(scenario,:) = [
            strength * 100, ...
            start_day, ...
            peak_cases, ...
            peak_day, ...
            cumulative_cases, ...
            peak_reduction, ...
            cumulative_reduction
        ];

        fprintf(['Intervention: %2.0f%% | Start day: %3d | ' ...
                 'Peak: %10.0f | Peak day: %3d | ' ...
                 'Cumulative: %12.0f\n'], ...
                 strength*100, start_day, ...
                 peak_cases, peak_day, cumulative_cases);
    end
end

% Save results to CSV
output_file = fullfile(project_dir, 'results', ...
                       'intervention_results.csv');

fid = fopen(output_file, 'w');

if fid == -1
    error('Could not create intervention results file.');
end

fprintf(fid, ...
    ['intervention_percent,start_day,peak_reported_cases,' ...
     'peak_day,cumulative_reported_cases,' ...
     'peak_reduction_percent,cumulative_reduction_percent\n']);

for i = 1:number_of_scenarios
    fprintf(fid, ...
        '%.2f,%d,%.4f,%d,%.4f,%.4f,%.4f\n', ...
        results(i,:));
end

fclose(fid);

% Plot intervention comparisons
figure;

subplot(2,1,1);
hold on;

for s = 1:length(intervention_strengths)
    strength = intervention_strengths(s) * 100;

    selected = results(:,1) == strength;

    plot(results(selected,2), results(selected,3), ...
         '-o', 'LineWidth', 1.5);
end

xlabel('Intervention start day');
ylabel('Peak daily reported cases');
title('Effect of Intervention on Peak Cases');
legend('0%', '20%', '40%', '60%', 'Location', 'best');
grid on;

subplot(2,1,2);
hold on;

for s = 1:length(intervention_strengths)
    strength = intervention_strengths(s) * 100;

    selected = results(:,1) == strength;

    plot(results(selected,2), results(selected,7), ...
         '-o', 'LineWidth', 1.5);
end

xlabel('Intervention start day');
ylabel('Cumulative case reduction (%)');
title('Effect of Intervention on Cumulative Cases');
legend('0%', '20%', '40%', '60%', 'Location', 'best');
grid on;

% Save figure
plot_file = fullfile(project_dir, 'results', ...
                     'intervention_analysis.png');

print(plot_file, '-dpng');

fprintf('\n------------------------------------\n');
fprintf('INTERVENTION ANALYSIS COMPLETED\n');
fprintf('------------------------------------\n');
fprintf('Scenarios simulated: %d\n', number_of_scenarios);
fprintf('Results saved: %s\n', output_file);
fprintf('Plot saved: %s\n', plot_file);
