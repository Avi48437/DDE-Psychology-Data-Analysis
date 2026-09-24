function output = summarize_matched_key_items_IPIPFMM(B1_matched, top_N)

% B1_matched:
% J x (K1 + 1) x number_of_runs
%
% The first column is the intercept.
% The remaining columns have already been matched across runs.

if nargin < 2
    top_N = 3;
end

survey_vars = [
    compose("EXT%d", 1:10), ...
    compose("EST%d", 1:10), ...
    compose("AGR%d", 1:10), ...
    compose("CSN%d", 1:10), ...
    compose("OPN%d", 1:10)];

[J, K1_plus_intercept, n_runs] = size(B1_matched);

K1 = K1_plus_intercept - 1;

all_key_items = strings(n_runs, K1, top_N);
all_key_strengths = zeros(n_runs, K1, top_N);

%% Extract key items from every matched run

for s = 1:n_runs

    B1 = B1_matched(:, :, s);

    active_cols = find(any(B1(:, 2:end) ~= 0, 1)) + 1;

    for k = 1:K1

        k_col = k + 1;
        other_cols = setdiff(active_cols, k_col);

        tmp = zeros(J, 1);

        for j = 1:J

            if isempty(other_cols)

                tmp(j) = max(0, B1(j, k_col));

            else

                loading_difference = ...
                    B1(j, k_col) - B1(j, other_cols);

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

            all_key_strengths(s, k, 1:n_keep) = ...
                sorted_tmp(1:n_keep);

        end

    end

end

%% Consensus top-N set for each matched factor

factor_names = "A1_" + string(1:K1);

consensus_items = strings(K1, top_N);
consensus_count = zeros(K1, 1);
consensus_percent = zeros(K1, 1);

for k = 1:K1

    factor_items = squeeze(all_key_items(:, k, :));

    if top_N == 1
        factor_items = factor_items(:);
    end

    % Ignore ordering within the top-N set
    sorted_items = sort(factor_items, 2);
    signatures = join(sorted_items, "|", 2);

    [unique_signatures, ~, group_index] = unique(signatures);
    signature_counts = accumarray(group_index, 1);

    [consensus_count(k), max_index] = max(signature_counts);

    consensus_percent(k) = ...
        100 * consensus_count(k) / n_runs;

    consensus_items(k, :) = ...
        split(unique_signatures(max_index), "|")';

end

ConsensusTable = array2table( ...
    consensus_items', ...
    'VariableNames', cellstr(factor_names), ...
    'RowNames', cellstr("Key Item " + string(1:top_N)));

%% Item-level appearance frequencies

factor_column = strings(0, 1);
item_column = strings(0, 1);
count_column = zeros(0, 1);
percent_column = zeros(0, 1);

for k = 1:K1

    factor_items = squeeze(all_key_items(:, k, :));
    factor_items = factor_items(:);
    factor_items = factor_items(strlength(factor_items) > 0);

    unique_items = unique(factor_items);

    for u = 1:numel(unique_items)

        item = unique_items(u);

        % Number of runs in which the item appears among top-N
        item_in_run = false(n_runs, 1);

        for s = 1:n_runs
            item_in_run(s) = any(all_key_items(s, k, :) == item);
        end

        item_count = sum(item_in_run);

        factor_column(end + 1, 1) = factor_names(k);
        item_column(end + 1, 1) = item;
        count_column(end + 1, 1) = item_count;
        percent_column(end + 1, 1) = 100 * item_count / n_runs;

    end

end

FrequencyTable = table( ...
    factor_column, ...
    item_column, ...
    count_column, ...
    percent_column, ...
    'VariableNames', ...
    {'Factor', 'Item', 'RunCount', 'RunPercent'});

FrequencyTable = sortrows( ...
    FrequencyTable, ...
    {'Factor', 'RunPercent'}, ...
    {'ascend', 'descend'});

%% Factor-level consensus summary

ConsensusSummary = table( ...
    factor_names', ...
    consensus_count, ...
    consensus_percent, ...
    'VariableNames', ...
    {'Factor', 'ConsensusSetCount', 'ConsensusSetPercent'});

%% Output

output.AllKeyItems = all_key_items;
output.AllKeyStrengths = all_key_strengths;
output.ConsensusTable = ConsensusTable;
output.ConsensusSummary = ConsensusSummary;
output.FrequencyTable = FrequencyTable;

end