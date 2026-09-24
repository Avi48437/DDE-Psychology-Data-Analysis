clear
clc
close all

%% Add folders
root_dir = "/Users/avinandanroy/Desktop/Research/DDE/DDE_Fitting_Pshycology2";

addpath(fullfile(root_dir, 'Scripts'));
addpath(fullfile(root_dir, 'Algorithms'));
addpath(fullfile(root_dir, 'Utilities'));
addpath(fullfile(root_dir, '2. IPIP 100', 'Fitting', 'Scripts'));

%% Load the IPIP-100 data
data_file = fullfile(root_dir, '2. IPIP 100', 'Data', 'B5.csv');

B5 = readtable( ...
    data_file, ...
    'FileType', 'text', ...
    'Delimiter', ',', ...
    'VariableNamingRule', 'preserve' ...
);

size(B5)
%% Encoding X
var_names = B5.Properties.VariableNames;

% Convert response labels to numeric values
X_text = strtrim(string(B5{:,:}));

[N, J] = size(X_text);
X_numeric = NaN(N, J);

X_numeric(X_text == "Very inaccurate") = 1;
X_numeric(X_text == "Moderately inaccurate") = 2;
X_numeric(X_text == "Neither inaccurate nor accurate") = 3;
X_numeric(X_text == "Moderately accurate") = 4;
X_numeric(X_text == "Very accurate") = 5;

% Store the numeric responses as a table
X = array2table( ...
    X_numeric, ...
    'VariableNames', var_names ...
);
%% Dimension
[N, J] = size(X);
R = NaN(N, J);

K1_max = floor(J/3); K2_max = floor(K1_max/3); epsilon = .00001;

%% ----------------------------------------------------
% 1) Convert each column to integer ranks
%    Missing values stay as NaN
% ----------------------------------------------------
for j = 1:J
    yj = X{:, j};   % extract raw column data

    % Convert nonnumeric columns to numeric codes
    if iscell(yj) || isstring(yj) || ischar(yj) || iscategorical(yj)
        yj = double(categorical(yj));
    end

    % Force to column vector
    yj = yj(:);

    % Check that we now have numeric data
    if ~isnumeric(yj)
        error('Column %d (%s) could not be converted to numeric.', j, var_names{j});
    end

    % Identify missing values
    isn = isnan(yj);

    % Rank distinct observed values in sorted order
    yj_nonan = yj(~isn);
    [~, ~, ic] = unique(yj_nonan, 'sorted');

    % Store ranks
    R(~isn, j) = ic;
end

%% ----------------------------------------------------
% 2) Number of rank levels per variable
% ----------------------------------------------------
Rlevels = zeros(1, J);
for j = 1:J
    if all(isnan(R(:, j)))
        Rlevels(j) = 0;
    else
        Rlevels(j) = max(R(:, j));
    end
end



X = X{:,:};
X(isnan(X)) = 3;

Z_init = simulate_gaussian_mixture(X,J);
[prop_in, B1_in, B2_in, gamma_in, A1_in, A2_in] = Normal_init_v2(Z_init, K1_max, K1_max, K2_max,epsilon);
%% Fit DDE once using the best hyperparameters

tau = 0.20;
t1 = 0.02;
t2 = 0.06;
temp = 0.70;

K1_max = floor(J / 3);
K2_max = floor(K1_max / 3);
epsilon = 0.00001;


fit_timer = tic;

[num_act1, num_act2, theta, theta2, ...
    gamma_CSP, p, q, A1_sample, A2_sample, ...
    Z, prop_CSP, B1_CSP, B2_CSP, ...
    A1_new, A2_new, it_num] = get_SAEM_RL_CSP( ...
    X, Z_init, R, Rlevels, ...
    prop_in, gamma_in, B1_in, B2_in, ...
    A1_in, A2_in, 1, 50, ...
    t1, t2, temp, tau);

fit_runtime = toc(fit_timer);

% Store fitted results

fit_results = struct();

fit_results.tau = tau;
fit_results.t1 = t1;
fit_results.t2 = t2;
fit_results.temp = temp;

fit_results.K1_max = K1_max;
fit_results.K2_max = K2_max;
fit_results.epsilon = epsilon;

fit_results.num_act1 = num_act1;
fit_results.num_act2 = num_act2;
fit_results.it_num = it_num;
fit_results.runtime = fit_runtime;

fit_results.theta = theta;
fit_results.theta2 = theta2;
fit_results.gamma_CSP = gamma_CSP;
fit_results.p = p;
fit_results.q = q;

fit_results.A1_sample = A1_sample;
fit_results.A2_sample = A2_sample;
fit_results.A1_new = A1_new;
fit_results.A2_new = A2_new;

fit_results.Z = Z;
fit_results.prop_CSP = prop_CSP;
fit_results.B1_CSP = B1_CSP;
fit_results.B2_CSP = B2_CSP;

fit_results.Z_init = Z_init;
fit_results.prop_in = prop_in;
fit_results.gamma_in = gamma_in;
fit_results.B1_in = B1_in;
fit_results.B2_in = B2_in;
fit_results.A1_in = A1_in;
fit_results.A2_in = A2_in;

fit_results.var_names = var_names;
% Display summary

fprintf('\nDDE fitting completed.\n');
fprintf('Active B1 dimensions: %d\n', num_act1);
fprintf('Active B2 dimensions: %d\n', num_act2);
fprintf('Number of iterations: %d\n', it_num);
fprintf('Runtime: %.2f minutes\n', fit_runtime / 60);
%fprintf('Results saved to:\n%s\n', results_file);

%% Save everything in one file

results = fit_results;

results_file = fullfile(root_dir, '2. IPIP 100', 'Fitting', 'Results', 'IPIP100_DDE_fit.mat');

save( ...
    results_file, ...
    'results', ...
    '-v7.3');

%% Plot fitted loading matrices

plot_loadings( ...
    fit_results.B1_CSP, ...
    fit_results.B2_CSP);