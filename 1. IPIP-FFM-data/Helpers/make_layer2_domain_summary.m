function summary_table = make_layer2_domain_summary(factor_cluster_table)

if isempty(factor_cluster_table)
    summary_table = table;
    return
end

A2_factors = unique( ...
    string(factor_cluster_table.DisplayA2), ...
    'stable');

summary_table = table;

%% ============================================================
% Construct summary table
%% ============================================================

for h = 1:numel(A2_factors)

    current_A2 = A2_factors(h);

    %% Positive side

    positive_rows = ...
        string(factor_cluster_table.DisplayA2) == current_A2 & ...
        string(factor_cluster_table.Direction) == "Positive";

    positive_table = factor_cluster_table(positive_rows,:);

    [positive_domains, ...
     positive_factors, ...
     positive_percent] = ...
        summarize_direction(positive_table);

    %% Negative side

    negative_rows = ...
        string(factor_cluster_table.DisplayA2) == current_A2 & ...
        string(factor_cluster_table.Direction) == "Negative";

    negative_table = factor_cluster_table(negative_rows,:);

    [negative_domains, ...
     negative_factors, ...
     negative_percent] = ...
        summarize_direction(negative_table);

    %% Number of rows needed for this Layer-2 column

    n_rows = max( ...
        numel(positive_domains), ...
        numel(negative_domains));

    if n_rows == 0
        continue
    end

    %% Pad positive side

    positive_domains(end+1:n_rows,1) = "";
    positive_factors(end+1:n_rows,1) = "";
    positive_percent(end+1:n_rows,1) = NaN;

    %% Pad negative side

    negative_domains(end+1:n_rows,1) = "";
    negative_factors(end+1:n_rows,1) = "";
    negative_percent(end+1:n_rows,1) = NaN;

    %% Add percentages to domain labels

    positive_labels = strings(n_rows,1);
    negative_labels = strings(n_rows,1);

    for j = 1:n_rows

        if strlength(positive_domains(j)) > 0

            positive_labels(j) = compose( ...
                "%s (%.1f%%)", ...
                positive_domains(j), ...
                positive_percent(j));
        end

        if strlength(negative_domains(j)) > 0

            negative_labels(j) = compose( ...
                "%s (%.1f%%)", ...
                negative_domains(j), ...
                negative_percent(j));
        end
    end

    active_column = repmat(current_A2,n_rows,1);

    current_summary = table( ...
        active_column, ...
        positive_labels, ...
        positive_factors, ...
        negative_labels, ...
        negative_factors, ...
        'VariableNames',{ ...
        'ActiveColumn', ...
        'PositiveDomain', ...
        'PositiveFactorsAssociated', ...
        'NegativeDomain', ...
        'NegativeFactorsAssociated'});

    summary_table = [summary_table; current_summary];
end


%% ============================================================
% Print formatted table without quotation marks
%% ============================================================

fprintf('\n');

fprintf('%-15s %-20s %-32s %-20s %-32s\n', ...
    'ActiveColumn', ...
    'PositiveDomain', ...
    'PositiveFactorsAssociated', ...
    'NegativeDomain', ...
    'NegativeFactorsAssociated');

separator = repmat('-',1,125);
fprintf('%s\n',separator);

previous_A2 = "";

for i = 1:height(summary_table)

    current_A2 = string(summary_table.ActiveColumn(i));

    %% Separator between different Layer-2 columns

    if i > 1 && current_A2 ~= previous_A2
        fprintf('%s\n',separator);
    end

    %% Print Layer-2 name only in the first row of its block

    if i == 1 || current_A2 ~= previous_A2
        active_label = current_A2;
    else
        active_label = "";
    end

    fprintf('%-15s %-20s %-32s %-20s %-32s\n', ...
        char(active_label), ...
        char(summary_table.PositiveDomain(i)), ...
        char(summary_table.PositiveFactorsAssociated(i)), ...
        char(summary_table.NegativeDomain(i)), ...
        char(summary_table.NegativeFactorsAssociated(i)));

    previous_A2 = current_A2;
end

fprintf('%s\n',separator);

end


function [domains,factors,percentages] = ...
    summarize_direction(direction_table)

if isempty(direction_table)

    domains = strings(0,1);
    factors = strings(0,1);
    percentages = zeros(0,1);

    return
end

domains = unique( ...
    string(direction_table.Domain), ...
    'stable');

n_domains = numel(domains);

factors = strings(n_domains,1);
percentages = zeros(n_domains,1);

for d = 1:n_domains

    current_domain = domains(d);

    idx = ...
        string(direction_table.Domain) == current_domain;

    domain_table = direction_table(idx,:);

    %% Average selection percentage for this domain

    percentages(d) = mean( ...
        domain_table.SelectionPercent);

    %% Sort associated factors by absolute B2 loading

    [~,order] = sort( ...
        abs(domain_table.MeanLoading), ...
        'descend');

    domain_table = domain_table(order,:);

    factors(d) = strjoin( ...
        string(domain_table.DisplayFactor), ...
        ", ");
end

%% Sort domains by selection percentage

[percentages,order] = sort( ...
    percentages, ...
    'descend');

domains = domains(order);
factors = factors(order);

end