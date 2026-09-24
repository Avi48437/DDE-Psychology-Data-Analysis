clear
clc
close all

%% Add folders

root_dir = '/Users/avinandanroy/Desktop/Research/DDE/DDE_Fitting_Pshycology2';

ipip_dir = fullfile(root_dir,'2. IPIP 100');
fit_results_dir = fullfile(ipip_dir,'Fitting','Results');
figures_dir = fullfile(ipip_dir,'Analysis','Figures');

addpath(fullfile(root_dir,'Scripts'));
addpath(fullfile(root_dir,'Algorithms'));
addpath(fullfile(root_dir,'Utilities'));

%% Load fitted IPIP-100 DDE model

load(fullfile(fit_results_dir,'IPIP100_DDE_fit.mat'));

B1_CSP = results.B1_CSP;
B2_CSP = results.B2_CSP;
var_names = string(results.var_names);

%% Plot and save compact loading figure

fig = figure( ...
    'Units', 'inches', ...
    'Position', [1 1 10 4.8], ...
    'Color', 'w', ...
    'Renderer', 'painters');

plot_loadings(B1_CSP, B2_CSP);

set(gcf, ...
    'Units', 'inches', ...
    'Position', [1 1 10 4.8], ...
    'Color', 'w', ...
    'Renderer', 'painters');

exportgraphics( ...
    gcf, ...
    fullfile(figures_dir,'IPIP100_Loadings.pdf'), ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'white');

%% ====================================================
% B1 anchor items
% =====================================================

top_N_B1 = 5;

J = size(B1_CSP, 1);

% Find active non-intercept B1 columns

active_factor_index = ...
    find(any(B1_CSP(:, 2:end) ~= 0, 1));

active_cols_B1 = active_factor_index + 1;
K1_active = length(active_cols_B1);

% Storage

top_items_B1 = strings(K1_active, top_N_B1);
top_index_B1 = zeros(K1_active, top_N_B1);
top_tmp_values_B1 = zeros(K1_active, top_N_B1);

% Find B1 anchors

for kk = 1:K1_active

    k_col = active_cols_B1(kk);
    other_cols = setdiff(active_cols_B1, k_col);

    tmp = zeros(J, 1);

    for j = 1:J

        if isempty(other_cols)

            tmp(j) = max(0, B1_CSP(j, k_col));

        else

            loading_differences = ...
                B1_CSP(j, k_col) - B1_CSP(j, other_cols);

            tmp(j) = max(0, min(loading_differences));

        end

    end

    [sorted_tmp, sort_idx] = sort(tmp, 'descend');

    positive_anchor_index = find(sorted_tmp > 0);
    n_keep = min(top_N_B1, length(positive_anchor_index));

    if n_keep > 0

        selected_idx = sort_idx(1:n_keep);

        top_index_B1(kk, 1:n_keep) = selected_idx;

        top_tmp_values_B1(kk, 1:n_keep) = ...
            sorted_tmp(1:n_keep);

        top_items_B1(kk, 1:n_keep) = ...
            var_names(selected_idx);

    end

end

% Print B1 anchors

for kk = 1:K1_active

    fprintf('\nB1 Active Factor %d:\n', ...
        active_factor_index(kk));

    positive_anchor_count = ...
        sum(top_tmp_values_B1(kk, :) > 0);

    if positive_anchor_count == 0

        fprintf('  No positive anchor items found.\n');

    else

        for m = 1:positive_anchor_count

            fprintf('  %d. %s (tmp = %.4f)\n', ...
                m, ...
                top_items_B1(kk, m), ...
                top_tmp_values_B1(kk, m));

        end

    end

end

%% ============================================================
% Construct anchor table for the first-layer latent factors
%% ============================================================

top_N = 5;

J = size(B1_CSP, 1);

% Active first-layer factors (excluding intercept)
active_factor_index = find(any(B1_CSP(:, 2:end) ~= 0, 1));
active_cols_B1 = active_factor_index + 1;

K1_active = length(active_cols_B1);

% Storage
top_items_B1 = strings(K1_active, top_N);
top_tmp_values = zeros(K1_active, top_N);

for kk = 1:K1_active

    current_col = active_cols_B1(kk);
    competing_cols = setdiff(active_cols_B1, current_col);

    tmp = zeros(J, 1);

    for j = 1:J

        if isempty(competing_cols)

            tmp(j) = max(0, B1_CSP(j, current_col));

        else

            loading_difference = ...
                B1_CSP(j, current_col) - ...
                B1_CSP(j, competing_cols);

            tmp(j) = max(0, min(loading_difference));

        end

    end

    [sorted_tmp, sort_idx] = sort(tmp, 'descend');

    positive_idx = find(sorted_tmp > 0);
    n_keep = min(top_N, length(positive_idx));

    if n_keep > 0

        selected_idx = sort_idx(1:n_keep);

        top_items_B1(kk, 1:n_keep) = ...
            string(var_names(selected_idx));

        top_tmp_values(kk, 1:n_keep) = ...
            sorted_tmp(1:n_keep);

    end

end

% Keep only interpretable factors

keep_factor = any(top_tmp_values > 0, 2);

anchor_codes = top_items_B1(keep_factor, :);

factor_index = active_factor_index(keep_factor);

factor_names = "A1_" + string(factor_index);

row_names = "Key Item " + string(1:top_N);

anchor_table = array2table( ...
    anchor_codes', ...
    'VariableNames', cellstr(factor_names), ...
    'RowNames', cellstr(row_names));

disp(anchor_table)

%% ============================================================
% Replace anchor variable codes (q1--q100) by survey questions
%% ============================================================

anchor_question_table = anchor_table;

var_mapping.Variable = string(var_mapping.Variable);
var_mapping.Question = string(var_mapping.Question);

for c = 1:width(anchor_question_table)

    for r = 1:height(anchor_question_table)

        question_code = string(anchor_table{r, c});

        if strlength(question_code) == 0
            continue
        end

        idx = find( ...
            var_mapping.Variable == question_code, ...
            1);

        if ~isempty(idx)

            anchor_question_table{r, c} = ...
                var_mapping.Question(idx);

        end

    end

end

%% ============================================================
% Display anchor questions in blocks of three factors
%% ============================================================

block_size = 3;

num_factors = width(anchor_question_table);

for first_col = 1:block_size:num_factors

    last_col = min( ...
        first_col + block_size - 1, ...
        num_factors);

    disp(anchor_question_table(:, first_col:last_col))

    fprintf('\n');

end
%% ====================================================
% B2 anchor factors
% =====================================================

top_N_B2 = 3;

K1_rows = size(B2_CSP, 1);

% All B1 factors are eligible for B2 anchor analysis

eligible_rows_B2 = 1:K1_rows;

% Find active non-intercept B2 columns

active_factor_index_B2 = ...
    find(any(B2_CSP(:, 2:end) ~= 0, 1));

active_cols_B2 = active_factor_index_B2 + 1;
K2_active = length(active_cols_B2);

% Storage

top_items_B2 = strings(K2_active, top_N_B2);
top_index_B2 = zeros(K2_active, top_N_B2);
top_tmp_values_B2 = zeros(K2_active, top_N_B2);

% Find B2 anchor factors

for kk = 1:K2_active

    k_col = active_cols_B2(kk);
    other_cols = setdiff(active_cols_B2, k_col);

    tmp = zeros(K1_rows, 1);

    for j = eligible_rows_B2

        if isempty(other_cols)

            tmp(j) = max(0, B2_CSP(j, k_col));

        else

            loading_differences = ...
                B2_CSP(j, k_col) - B2_CSP(j, other_cols);

            tmp(j) = max(0, min(loading_differences));

        end

    end

    [sorted_tmp, local_sort_idx] = ...
        sort(tmp(eligible_rows_B2), 'descend');

    positive_anchor_index = find(sorted_tmp > 0);

    n_keep = min( ...
        top_N_B2, ...
        length(positive_anchor_index));

    if n_keep > 0

        selected_local_idx = ...
            local_sort_idx(1:n_keep);

        selected_rows = ...
            eligible_rows_B2(selected_local_idx);

        top_index_B2(kk, 1:n_keep) = ...
            selected_rows;

        top_tmp_values_B2(kk, 1:n_keep) = ...
            sorted_tmp(1:n_keep);

        top_items_B2(kk, 1:n_keep) = ...
            "A1_" + string(selected_rows);

    end

end

% Print B2 anchor factors

for kk = 1:K2_active

    fprintf('\nB2 Active Factor %d:\n', ...
        active_factor_index_B2(kk));

    positive_anchor_count = ...
        sum(top_tmp_values_B2(kk, :) > 0);

    if positive_anchor_count == 0

        fprintf('  No positive anchor factors found.\n');

    else

        for m = 1:positive_anchor_count

            fprintf('  %d. %s (tmp = %.4f)\n', ...
                m, ...
                top_items_B2(kk, m), ...
                top_tmp_values_B2(kk, m));

        end

    end

end

%% ====================================================
% Display B2 loading tables
% =====================================================

for kk = 1:K2_active

    factor_number = active_factor_index_B2(kk);
    loading_column = active_cols_B2(kk);

    B2_loading_table = table( ...
        (1:K1_rows)', ...
        B2_CSP(:, loading_column), ...
        'VariableNames', ...
        {'A1_Factor', 'B2_Loading'});

    B2_loading_table = sortrows( ...
        B2_loading_table, ...
        'B2_Loading', ...
        'descend');

    fprintf('\nB2 Active Factor %d Loadings:\n', ...
        factor_number);

    disp(B2_loading_table)

end