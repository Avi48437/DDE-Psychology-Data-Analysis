clear
clc
close all
%% Add folders
root_dir = '/Users/avinandanroy/Desktop/Research/DDE/DDE_Fitting_Pshycology2';

fmm_dir = fullfile(root_dir,'1. IPIP-FFM-data');
data_dir = fullfile(fmm_dir,'Data');
fit_results_dir = fullfile(fmm_dir,'Fitting','Results');
helper_dir = fullfile(fmm_dir,'Helpers');

addpath(fullfile(root_dir,'Scripts'));
addpath(fullfile(root_dir,'Algorithms'));
addpath(fullfile(root_dir,'Utilities'));
addpath(helper_dir);

load(fullfile(fit_results_dir,'IPIP_DDE_100_results6_comp.mat'))

%% Load the IPIP-FFM data
data_file = fullfile(data_dir,'data-final.csv');

B5 = readtable( ...
    data_file, ...
    'FileType', 'text', ...
    'Delimiter', '\t' ...
    );
% Load the cleaned IPIP-FFM data

reduced_data_file = fullfile(data_dir,'B5_Red.mat');

load(reduced_data_file, 'B5_Red');

fprintf('Original data size: %d rows x %d variables\n', ...
    height(B5), width(B5));

fprintf('Cleaned data size : %d rows x %d variables\n', ...
    height(B5_Red), width(B5_Red));
%% Collapse countries with fewer than 1000 observations into "OTHER"

B5CM = B5_Red;

country = categorical(B5_Red.country);

country_names = categories(country);
country_counts = countcats(country);

large_countries = country_names(country_counts >= 1000);

country_new = string(B5_Red.country);
country_new(~ismember(country_new, string(large_countries))) = "OTHER";

B5CM.country = categorical(country_new);


% Keep respondents with valid 1--5 responses on all 50 items
item_names = ["EXT"+string(1:10),...
    "EST"+string(1:10),...
    "AGR"+string(1:10),...
    "CSN"+string(1:10),...
    "OPN"+string(1:10)];


X_items = B5_Red{:,item_names};
%%
complete_rows = all(ismember(X_items,1:5),2);
B5_Red_Comp = B5_Red(complete_rows,:);
X_items = B5_Red{:,item_names};
complete_rows = all(ismember(X_items,1:5),2);
B5_Red_Comp = B5_Red(complete_rows,:);

% Collapse countries with fewer than 1000 respondents
B5CM = B5_Red_Comp;

country = categorical(B5_Red_Comp.country);
country_names = categories(country);
country_counts = countcats(country);
large_countries = country_names(country_counts >= 1000);

country_new = string(B5_Red_Comp.country);
country_new(~ismember(country_new,string(large_countries))) = "OTHER";
B5CM.country = categorical(country_new);

fprintf('Complete-case dataset: %d respondents\n',height(B5CM));
%% ============================================================
% Post-fit analysis
% ============================================================

%% 1. Completed repetitions

valid_iterations = find(completed);

if isempty(valid_iterations)
    error('No completed repetitions were found.');
end

fprintf('Completed repetitions: %d\n', length(valid_iterations));

%% ============================================================
% 2. Frequency of active-dimension combinations
% ============================================================

n_valid = length(valid_iterations);

K1_active = zeros(n_valid, 1);
K2_active = zeros(n_valid, 1);

for i = 1:n_valid

    r = valid_iterations(i);

    B1_r = results(r).B1_CSP;
    B2_r = results(r).B2_CSP;

    % Count active non-intercept loading columns
    K1_active(i) = sum(any(B1_r(:, 2:end) ~= 0, 1));
    K2_active(i) = sum(any(B2_r(:, 2:end) ~= 0, 1));

end

dimension_pairs = [K1_active, K2_active];

[unique_pairs, ~, group_id] = unique( dimension_pairs,'rows');

frequency = accumarray(group_id, 1);
percent = 100 * frequency / n_valid;

dimension_summary = table( ...
    unique_pairs(:, 1), ...
    unique_pairs(:, 2), ...
    frequency, ...
    percent, ...
    'VariableNames', ...
    {'K1_Active', 'K2_Active', 'Count', 'Percent'});

dimension_summary = sortrows( ...
    dimension_summary, ...
    {'Count', 'K1_Active', 'K2_Active'}, ...
    {'descend', 'ascend', 'ascend'});

disp('Active-dimension combinations:')
disp(dimension_summary)

%% Summary of active dimensions

active_dimension_summary = table( ...
    [mean(K1_active); mean(K2_active)], ...
    [std(K1_active); std(K2_active)], ...
    [median(K1_active); median(K2_active)], ...
    [min(K1_active); min(K2_active)], ...
    [max(K1_active); max(K2_active)], ...
    'VariableNames', ...
    {'Mean', 'SD', 'Median', 'Minimum', 'Maximum'}, ...
    'RowNames', ...
    {'Layer 1', 'Layer 2'});

disp('Summary of active dimensions:')
disp(active_dimension_summary)

%% ============================================================
% 3. Average matched loading matrices
% ============================================================

average_fit = average_matched_loadings(results, 0.1);

average_B1 = average_fit.average_B1;
average_B2 = average_fit.average_B2;

plot_loadings(average_B1, average_B2);

%% ============================================================
% 4. Extract IPIP-FFM key items from average B1
% ============================================================

top_N_B1 = 3;

[anchor_table, ...
 anchor_question_table, ...
 anchor_strength] = ...
    find_key_item_IPIPFMM(average_B1, top_N_B1);

disp('Key-item variable names:')
disp(anchor_table)

disp('Key-item question text:')
disp(anchor_question_table)

%% Display question table in blocks

block_size = 3;
num_factors = width(anchor_question_table);

for start_col = 1:block_size:num_factors

    end_col = min( ...
        start_col + block_size - 1, ...
        num_factors);

    disp(anchor_question_table(:, start_col:end_col))
    fprintf('\n');

end

%% ============================================================
% 5. Inspect active second-layer loadings
% ============================================================

active_B2_columns = find( ...
    any(average_B2(:, 2:end) ~= 0, 1));

if isempty(active_B2_columns)

    fprintf('The average fit has no active second-layer factors.\n');

    B2_loading_table = table();

else

    A1_factor_index = (1:size(average_B2, 1))';

    B2_loading_table = array2table( ...
        average_B2(:, active_B2_columns + 1), ...
        'VariableNames', ...
        cellstr("A2_" + string(active_B2_columns)));

    B2_loading_table = addvars( ...
        B2_loading_table, ...
        A1_factor_index, ...
        'Before', 1, ...
        'NewVariableNames', 'A1_Factor');

    disp('Active average B2 loadings:')
    disp(B2_loading_table)

end
%% ============================================================
% 6. Assign IDs and create matched groups
% ============================================================

n_combinations = height(dimension_summary);

dimension_summary.CombinationID = ...
    "C" + string((1:n_combinations)');

dimension_summary = movevars( ...
    dimension_summary, ...
    'CombinationID', ...
    'Before', 1);

disp('Active-dimension combinations with IDs:')
disp(dimension_summary)

for c = 1:n_combinations

    id = dimension_summary.CombinationID(c);

    target_K1 = dimension_summary.K1_Active(c);
    target_K2 = dimension_summary.K2_Active(c);

    repetition_ids = valid_iterations( ...
        K1_active == target_K1 & ...
        K2_active == target_K2);

    group_results = results(repetition_ids);

    matched_group = match_loading_group( ...
        group_results, repetition_ids);

    eval(id + " = matched_group;");

end
%%
combination_id = "C1";

C = eval(char(combination_id));

C = print_IPIPFMM_combination( ...
    C, ...
    combination_id, ...
    3, ...       % top key items
    1e-10, ...   % zero tolerance
    12);          % number of Layer-1 factors to display and use in Layer 2

assignin('base',char(combination_id),C);

%% Plot averaged matched loadings for one combination

combination_id = "C1";

row_idx = find(dimension_summary.CombinationID == combination_id, 1);
target_K1 = dimension_summary.K1_Active(row_idx);
target_K2 = dimension_summary.K2_Active(row_idx);

repetition_ids = valid_iterations(K1_active == target_K1 & K2_active == target_K2);

average_fit_comb = average_matched_loadings(results(repetition_ids), 0.1);

average_B1_comb = average_fit_comb.average_B1;
average_B2_comb = average_fit_comb.average_B2;

plot_loadings(average_B1_comb, average_B2_comb);
%%
top_N = 3;

[anchor_table_comb, ...
    anchor_question_table_comb, ...
    anchor_strength_comb] = ...
    find_key_item_IPIPFMM(average_B1_comb,top_N);

% Same factor ordering as B1 plot

groups_FMM = {[1 14],[2 12],[4 15],[6 3],[11 7 5 9]};
trait_labels = ["EXT","EST","AGR","CSN","OPN"];

factor_order = [groups_FMM{:}];

factor_names = "A1_" + string(factor_order);

% Find these columns in the anchor table

current_names = string(anchor_question_table_comb.Properties.VariableNames);

[found,order_idx] = ismember(factor_names,current_names);

if any(~found)
    error('Could not find factors: %s',strjoin(factor_names(~found),', '));
end

anchor_question_table_print = anchor_question_table_comb(:,order_idx);

% Domain labels in the same group order

domains = strings(1,numel(factor_order));

pos = 1;

for d = 1:numel(groups_FMM)
    n_d = numel(groups_FMM{d});
    domains(pos:pos+n_d-1) = trait_labels(d);
    pos = pos+n_d;
end

% Rename columns

display_names = factor_names + " (" + domains + ")";

anchor_question_table_print.Properties.VariableNames = cellstr(display_names);

% Print 3 factors at a time = 4 blocks

factors_per_block = 3;
n_factors = width(anchor_question_table_print);
n_blocks = ceil(n_factors/factors_per_block);

for b = 1:n_blocks

    first_col = (b-1)*factors_per_block + 1;
    last_col = min(b*factors_per_block,n_factors);

    fprintf('\n============================================================\n');
    fprintf('Layer-1 factors: block %d of %d\n',b,n_blocks);
    fprintf('============================================================\n');

    disp(anchor_question_table_print(:,first_col:last_col))
end
%% ============================================================
% Compare IPIP-FMM C1 (12,1) and IPIP-100 (8,1)
% B1 and B2 in one figure
% ============================================================
groups_FMM = {[1 14],[2 12],[4 15],[6 3],[11 7 5 9]};
groups_100 = {[1 2],[8 18],[13],[5 32],[10]};

block_sizes_FMM = [10 10 10 10 10];
block_sizes_100 = [20 20 20 20 20];

trait_labels = {'EXT','EST','AGR','CSN','OPN'};
plot_loadings_comparison2( ...
    average_B1_comb,average_B2_comb, ...
    groups_FMM,block_sizes_FMM,'IPIP-FMM C1 (12,1)', ...
    B1_IPIP100_grouped,B2_IPIP100, ...
    groups_100,block_sizes_100,'IPIP-100 (8,1)', ...
    trait_labels);