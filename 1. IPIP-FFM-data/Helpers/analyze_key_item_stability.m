function output = analyze_key_item_stability( ...
    results, valid_iterations, K1_active, K2_active, ...
    target_K1, target_K2, top_N)

%% Identify repetitions

repetition_ids = valid_iterations( ...
    K1_active == target_K1 & K2_active == target_K2);

n_runs = numel(repetition_ids);

if n_runs == 0
    error('No runs found for K1 = %d and K2 = %d.', ...
        target_K1, target_K2);
end

%% Obtain matched average and use it as the reference

results_subset = results(repetition_ids);
average_fit = average_matched_loadings(results_subset, 0);
reference_B1 = average_fit.average_B1;

reference_active = find(any(reference_B1(:, 2:end) ~= 0, 1)) + 1;
n_factors = numel(reference_active);

if n_factors ~= target_K1
    warning('Reference average contains %d active factors, expected %d.', ...
        n_factors, target_K1);
end

%% IPIP variable names

survey_vars = [
    compose("EXT%d", 1:10), ...
    compose("EST%d", 1:10), ...
    compose("AGR%d", 1:10), ...
    compose("CSN%d", 1:10), ...
    compose("OPN%d", 1:10)];

%% Storage

all_key_items = strings(n_runs, n_factors, top_N);
matched_B1 = cell(n_runs, 1);

%% Match every run to the average reference

for s = 1:n_runs

    B1 = results(repetition_ids(s)).B1_CSP;
    run_active = find(any(B1(:, 2:end) ~= 0, 1)) + 1;

    if numel(run_active) ~= n_factors
        error('Repetition %d has %d active factors, expected %d.', ...
            repetition_ids(s), numel(run_active), n_factors);
    end

    reference_loadings = reference_B1(:, reference_active);
    run_loadings = B1(:, run_active);

    similarity = zeros(n_factors);

    for a = 1:n_factors
        for b = 1:n_factors

            x = reference_loadings(:, a);
            y = run_loadings(:, b);

            similarity(a, b) = ...
                dot(x, y) / (norm(x) * norm(y) + eps);
        end
    end

    assignment = matchpairs(-similarity, 1e10);

    aligned_B1 = zeros(size(B1));
    aligned_B1(:, 1) = B1(:, 1);

    for a = 1:size(assignment, 1)

        reference_position = assignment(a, 1);
        run_position = assignment(a, 2);

        aligned_B1(:, reference_active(reference_position)) = ...
            B1(:, run_active(run_position));
    end

    matched_B1{s} = aligned_B1;

    %% Extract key items from each aligned factor

    for k = 1:n_factors

        k_col = reference_active(k);
        other_cols = setdiff(reference_active, k_col);

        tmp = zeros(size(B1, 1), 1);

        for j = 1:size(B1, 1)

            if isempty(other_cols)
                tmp(j) = max(0, aligned_B1(j, k_col));
            else
                loading_difference = ...
                    aligned_B1(j, k_col) - aligned_B1(j, other_cols);

                tmp(j) = max(0, min(loading_difference));
            end
        end

        [sorted_tmp, sort_index] = sort(tmp, 'descend');
        positive_index = find(sorted_tmp > 0);
        n_keep = min(top_N, numel(positive_index));

        if n_keep > 0
            selected_items = sort_index(1:n_keep);

            all_key_items(s, k, 1:n_keep) = ...
                survey_vars(selected_items);
        end
    end
end

%% Compute exact agreement across repetitions

factor_name = "A1_" + string(reference_active - 1);

top1_item = strings(n_factors, 1);
top1_frequency = zeros(n_factors, 1);
top1_percent = zeros(n_factors, 1);

exact_order_frequency = zeros(n_factors, 1);
exact_order_percent = zeros(n_factors, 1);

exact_set_frequency = zeros(n_factors, 1);
exact_set_percent = zeros(n_factors, 1);

consensus_items = strings(n_factors, top_N);

for k = 1:n_factors

    factor_items = squeeze(all_key_items(:, k, :));

    if top_N == 1
        factor_items = factor_items(:);
    end

    %% Most frequent top item

    top1 = factor_items(:, 1);
    [unique_top1, ~, top1_group] = unique(top1);
    top1_counts = accumarray(top1_group, 1);

    [top1_frequency(k), max_index] = max(top1_counts);
    top1_item(k) = unique_top1(max_index);
    top1_percent(k) = 100 * top1_frequency(k) / n_runs;

    %% Most frequent ordered top-N combination

    ordered_signature = join(factor_items, "|", 2);
    [unique_ordered, ~, ordered_group] = unique(ordered_signature);
    ordered_counts = accumarray(ordered_group, 1);

    [exact_order_frequency(k), max_index] = max(ordered_counts);
    exact_order_percent(k) = ...
        100 * exact_order_frequency(k) / n_runs;

    consensus_items(k, :) = ...
        split(unique_ordered(max_index), "|")';

    %% Most frequent unordered top-N set

    sorted_items = sort(factor_items, 2);
    set_signature = join(sorted_items, "|", 2);

    [unique_sets, ~, set_group] = unique(set_signature);
    set_counts = accumarray(set_group, 1);

    exact_set_frequency(k) = max(set_counts);
    exact_set_percent(k) = ...
        100 * exact_set_frequency(k) / n_runs;
end

%% Summary table

SummaryTable = table( ...
    factor_name', ...
    top1_item, ...
    top1_frequency, ...
    top1_percent, ...
    exact_order_frequency, ...
    exact_order_percent, ...
    exact_set_frequency, ...
    exact_set_percent, ...
    'VariableNames', { ...
    'Factor', ...
    'MostFrequentTopItem', ...
    'Top1Count', ...
    'Top1Percent', ...
    'ExactOrderedTopNCount', ...
    'ExactOrderedTopNPercent', ...
    'ExactTopNSetCount', ...
    'ExactTopNSetPercent'});

%% Consensus key-item table

ConsensusTable = array2table( ...
    consensus_items', ...
    'VariableNames', cellstr(factor_name), ...
    'RowNames', cellstr("Key Item " + string(1:top_N)));

%% Store output

output.RepetitionIDs = repetition_ids;
output.ReferenceB1 = reference_B1;
output.MatchedB1 = matched_B1;
output.AllKeyItems = all_key_items;
output.SummaryTable = SummaryTable;
output.ConsensusTable = ConsensusTable;

fprintf('\nK1 = %d, K2 = %d, repetitions = %d\n', ...
    target_K1, target_K2, n_runs);

disp('Consensus key-item table:')
disp(ConsensusTable)

end