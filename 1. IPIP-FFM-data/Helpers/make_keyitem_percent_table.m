function [keyitem_table, factor_mapping] = ...
    make_keyitem_percent_table(stability_output, top_N)

if nargin < 2
    top_N = 3;
end

F = stability_output.FrequencyTable;

factors = unique(F.Factor, 'stable');
n_factors = numel(factors);

display_values = strings(top_N, n_factors);
factor_domain = strings(n_factors, 1);

for k = 1:n_factors

    rows = F(F.Factor == factors(k), :);
    rows = sortrows(rows, 'RunPercent', 'descend');

    n_keep = min(top_N, height(rows));

    for q = 1:n_keep

        display_values(q, k) = ...
            rows.Item(q) + ...
            " (" + ...
            compose("%.1f", rows.RunPercent(q)) + ...
            "%)";

    end

    % Determine domain from the most frequent key item
    if n_keep > 0
        factor_domain(k) = extractBefore(rows.Item(1), 4);
    else
        factor_domain(k) = "OTHER";
    end

end

%% Order factors by personality domain

domain_order = ["EXT", "OPN", "AGR", "CSN", "EST", "OTHER"];
domain_rank = zeros(n_factors, 1);

for k = 1:n_factors

    idx = find(domain_order == factor_domain(k), 1);

    if isempty(idx)
        domain_rank(k) = numel(domain_order);
    else
        domain_rank(k) = idx;
    end

end

% Preserve the original factor order within each domain
original_order = (1:n_factors)';

[~, factor_order] = sortrows( ...
    [domain_rank, original_order], ...
    [1, 2]);

display_values = display_values(:, factor_order);

%% Rename columns as Factor1, ..., FactorK

new_factor_names = "Factor" + string(1:n_factors);
row_names = "Item" + string(1:top_N);

keyitem_table = array2table( ...
    display_values, ...
    'VariableNames', cellstr(new_factor_names), ...
    'RowNames', cellstr(row_names));

%% Mapping from display factors to original matched factors

factor_mapping = table( ...
    new_factor_names', ...
    factors(factor_order), ...
    factor_domain(factor_order), ...
    'VariableNames', ...
    {'DisplayFactor', 'OriginalFactor', 'Domain'});

end