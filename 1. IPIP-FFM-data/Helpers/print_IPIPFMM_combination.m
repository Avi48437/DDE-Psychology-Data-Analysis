function C = print_IPIPFMM_combination( ...
    C,combination_id,top_N,zero_tol,n_factors_to_show)

%% ============================================================
% Defaults
%% ============================================================

if nargin < 2 || isempty(combination_id)
    combination_id = "Combination";
end

if nargin < 3 || isempty(top_N)
    top_N = 3;
end

if nargin < 4 || isempty(zero_tol)
    zero_tol = 1e-10;
end

if nargin < 5 || isempty(n_factors_to_show)
    n_factors_to_show = Inf;
end

combination_id = string(combination_id);

factors_per_block = 3;

min_abs_loading = 0.05;
min_selection_percent = 50;
min_sign_consistency = 70;

required_fields = {'B1','B2','K1','K2'};

for j = 1:numel(required_fields)
    if ~isfield(C,required_fields{j})
        error('C is missing the field %s.',required_fields{j});
    end
end

if isfield(C,'RepetitionIDs')
    n_runs = numel(C.RepetitionIDs);
else
    n_runs = size(C.B1,3);
end

%% ============================================================
% Header
%% ============================================================

fprintf('\n========================================\n');
fprintf('%s: K1 = %d, K2 = %d, runs = %d\n', ...
    combination_id,C.K1,C.K2,n_runs);
fprintf('========================================\n');

%% ============================================================
% 1. Jointly align B1 and B2 across repetitions
%% ============================================================

[B1_aligned_full,B2_aligned_full,alignment_info] = ...
    align_combination_B1_B2(C.B1,C.B2,C.K1,C.K2,zero_tol);

C.B1AlignedFull = B1_aligned_full;
C.B2AlignedFull = B2_aligned_full;
C.AlignmentInfo = alignment_info;

%% ============================================================
% 2. Determine ORIGINAL factor identities from reference run
%
% align_combination_B1_B2 uses repetition 1 as the reference.
% The aligned factors 1:K1 are compressed positions, NOT the
% original A1 column numbers.
%% ============================================================

reference_iteration = 1;

B1_reference = C.B1(:,:,reference_iteration);
B2_reference = C.B2(:,:,reference_iteration);

reference_A1 = find(any(abs(B1_reference(:,2:end)) > zero_tol,1));
reference_A2 = find(any(abs(B2_reference(:,2:end)) > zero_tol,1));

if numel(reference_A1) ~= C.K1
    error('Reference run has %d active A1 factors but C.K1 = %d.', ...
        numel(reference_A1),C.K1);
end

if numel(reference_A2) ~= C.K2
    error('Reference run has %d active A2 factors but C.K2 = %d.', ...
        numel(reference_A2),C.K2);
end

C.ReferenceActiveA1 = reference_A1;
C.ReferenceActiveA2 = reference_A2;

fprintf('Original active A1 columns in reference run: [%s]\n',num2str(reference_A1));
fprintf('Original active A2 columns in reference run: [%s]\n',num2str(reference_A2));

%% ============================================================
% 3. Layer-1 key-item stability
%% ============================================================

C.KeyItemStability = ...
    summarize_matched_key_items_IPIPFMM(B1_aligned_full,top_N);

[key_item_table,factor_mapping] = ...
    make_keyitem_percent_table(C.KeyItemStability,top_N);

n_factors_total = width(key_item_table);
n_key_items = height(key_item_table);

if n_factors_total ~= C.K1
    error('Key-item table contains %d factors but C.K1 = %d.', ...
        n_factors_total,C.K1);
end

old_display_names = string(key_item_table.Properties.VariableNames);

%% ============================================================
% 4. Compute Layer-1 factor consistency
%% ============================================================

mean_consistency = zeros(n_factors_total,1);
minimum_consistency = zeros(n_factors_total,1);

for j = 1:n_factors_total

    labels = string(key_item_table{:,j});
    percentages = nan(n_key_items,1);

    for i = 1:n_key_items

        current_label = labels(i);

        if ismissing(current_label) || strlength(current_label)==0
            continue
        end

        token = regexp(char(current_label),'\(([\d.]+)%\)','tokens','once');

        if ~isempty(token)
            percentages(i) = str2double(token{1});
        end
    end

    percentages = percentages(~isnan(percentages));

    if isempty(percentages)
        mean_consistency(j) = 0;
        minimum_consistency(j) = 0;
    else
        mean_consistency(j) = mean(percentages);
        minimum_consistency(j) = min(percentages);
    end
end

ranking_table = table( ...
    (1:n_factors_total)', ...
    old_display_names', ...
    mean_consistency, ...
    minimum_consistency, ...
    'VariableNames',{ ...
    'AlignedColumn', ...
    'AlignedDisplayFactor', ...
    'MeanConsistency', ...
    'MinimumConsistency'});

ranking_table = sortrows( ...
    ranking_table, ...
    {'MeanConsistency','MinimumConsistency','AlignedColumn'}, ...
    {'descend','descend','ascend'});

%% ============================================================
% 5. Retain requested number of factors
%% ============================================================

if isinf(n_factors_to_show)

    n_factors_selected = n_factors_total;

else

    if ~isscalar(n_factors_to_show) || ...
            ~isfinite(n_factors_to_show) || ...
            n_factors_to_show < 1

        error('n_factors_to_show must be a positive integer or Inf.');
    end

    n_factors_selected = min(floor(n_factors_to_show),n_factors_total);
end

ranking_table = ranking_table(1:n_factors_selected,:);

% IMPORTANT:
% These are positions in the compressed aligned matrices.
aligned_factor_rows = ranking_table.AlignedColumn;

fprintf('\nDisplaying the %d most consistent of %d Layer-1 factors.\n', ...
    n_factors_selected,n_factors_total);

%% ============================================================
% 6. Reorder key-item table and factor mapping
%% ============================================================

ordered_old_names = old_display_names(aligned_factor_rows);

key_item_table = key_item_table(:,aligned_factor_rows);

mapping_names = string(factor_mapping.DisplayFactor);

[found,mapping_order] = ismember(ordered_old_names,mapping_names);

if any(~found)
    error('Could not map all aligned factors.');
end

factor_mapping = factor_mapping(mapping_order,:);

%% ============================================================
% 7. Map aligned slots back to ORIGINAL DDE column numbers
%
% Example:
%
% aligned slot 8 may correspond to original A1_11.
%
% This is the part that was wrong before.
%% ============================================================

original_A1_indices = reference_A1(aligned_factor_rows);
original_A1_names = "A1_" + string(original_A1_indices(:));

factor_mapping.DisplayFactor = original_A1_names;
factor_mapping.OriginalFactor = original_A1_names;

%% ============================================================
% 8. Restrict aligned matrices
%
% Use ALIGNED positions for indexing the compressed matrices.
% Use ORIGINAL names only for labels.
%% ============================================================

C.B1Aligned = B1_aligned_full( ...
    :, ...
    [1,aligned_factor_rows(:)'+1], ...
    :);

C.B2Aligned = B2_aligned_full( ...
    aligned_factor_rows, ...
    :, ...
    :);

C.SelectedAlignedFactorRows = aligned_factor_rows;
C.SelectedOriginalFactorIndices = original_A1_indices(:);
C.SelectedOriginalFactors = original_A1_names;

%% ============================================================
% 9. Store tables using true original names
%% ============================================================

key_item_table.Properties.VariableNames = cellstr(original_A1_names);

C.KeyItemTable = key_item_table;
C.FactorMapping = factor_mapping;

C.FactorConsistency = table( ...
    original_A1_names, ...
    aligned_factor_rows(:), ...
    string(factor_mapping.Domain), ...
    ranking_table.MeanConsistency, ...
    ranking_table.MinimumConsistency, ...
    'VariableNames',{ ...
    'OriginalFactor', ...
    'AlignedPosition', ...
    'Domain', ...
    'MeanConsistency', ...
    'MinimumConsistency'});

%% ============================================================
% 10. Question-text table
%% ============================================================

C.KeyItemQuestionTable = ...
    create_IPIP_question_table(C.KeyItemTable,C.FactorMapping);

%% ============================================================
% 11. Print Layer-1 factors
%% ============================================================

fprintf('\nLayer-1 key-item stability:\n');

n_blocks = ceil(n_factors_selected/factors_per_block);

for b = 1:n_blocks

    first_column = (b-1)*factors_per_block + 1;
    last_column = min(b*factors_per_block,n_factors_selected);

    fprintf('\n============================================================\n');
    fprintf('Layer-1 factors: block %d of %d\n',b,n_blocks);
    fprintf('============================================================\n');

    disp(C.KeyItemQuestionTable(:,first_column:last_column))
end

%% ============================================================
% 12. Average selected aligned matrices
%% ============================================================

C.AverageB1 = mean(C.B1Aligned,3);
C.AverageB2 = mean(C.B2Aligned,3);

%% ============================================================
% 13. No Layer-2 factors
%% ============================================================

if C.K2 == 0

    C.AverageB2Table = table;
    C.SecondLayerTable = table;
    C.Layer2FactorClusters = table;
    C.Layer2DomainSummary = table;

    fprintf('\n%s has no active Layer-2 factors.\n',combination_id);

    return
end

%% ============================================================
% 14. Layer-2 interpretation
%
% Pass ORIGINAL A2 identities as well.
%% ============================================================

[C.AverageB2Table,C.SecondLayerTable] = ...
    create_layer2_interpretation_tables( ...
        C.B2Aligned, ...
        C.FactorMapping, ...
        reference_A2, ...
        zero_tol);

fprintf('\nAverage Layer-2 loadings for selected Layer-1 factors:\n');
disp(C.AverageB2Table)

%% ============================================================
% 15. Stable Layer-2 relationships
%% ============================================================

keep = ...
    abs(C.SecondLayerTable.MeanLoading) >= min_abs_loading & ...
    C.SecondLayerTable.SelectionPercent >= min_selection_percent & ...
    C.SecondLayerTable.SignConsistencyPercent >= min_sign_consistency;

C.Layer2FactorClusters = C.SecondLayerTable(keep,:);

if isempty(C.Layer2FactorClusters)

    C.Layer2DomainSummary = table;

    fprintf('\nNo Layer-2 relationships passed the stability thresholds.\n');

    return
end

direction = repmat("Positive",height(C.Layer2FactorClusters),1);
direction(C.Layer2FactorClusters.MeanLoading < 0) = "Negative";

C.Layer2FactorClusters = addvars( ...
    C.Layer2FactorClusters, ...
    direction, ...
    'After','DisplayA2', ...
    'NewVariableNames','Direction');

C.Layer2FactorClusters.AbsMeanLoading = ...
    abs(C.Layer2FactorClusters.MeanLoading);

C.Layer2FactorClusters = sortrows( ...
    C.Layer2FactorClusters, ...
    {'DisplayA2','Direction','AbsMeanLoading'}, ...
    {'ascend','ascend','descend'});

%% ============================================================
% 16. Positive-versus-negative domain summary
%% ============================================================

C.Layer2DomainSummary = ...
    create_layer2_domain_summary(C.Layer2FactorClusters);

print_layer2_domain_summary(C.Layer2DomainSummary);

end

function question_table = ...
    create_IPIP_question_table(key_item_table,factor_mapping)

lookup = IPIP_item_lookup();

n_items = height(key_item_table);
n_factors = width(key_item_table);

question_matrix = strings(n_items,n_factors);

for j = 1:n_factors

    for i = 1:n_items

        current_label = string(key_item_table{i,j});

        if ismissing(current_label) || strlength(current_label)==0
            continue
        end

        token = regexp( ...
            char(current_label), ...
            '^([A-Z]{3}\d+)\s*\(([\d.]+)%\)$', ...
            'tokens','once');

        if isempty(token)
            question_matrix(i,j) = current_label;
            continue
        end

        item_code = string(token{1});
        stability_percent = string(token{2});

        item_position = find(lookup(:,1)==item_code,1);

        if isempty(item_position)

            question_matrix(i,j) = current_label;

        else

            question_text = lookup(item_position,2);

            question_matrix(i,j) = ...
                question_text + ...
                " [" + item_code + "] " + ...
                "(" + stability_percent + "%)";
        end
    end
end

column_names = strings(1,n_factors);

for j = 1:n_factors
    factor_name = string(factor_mapping.DisplayFactor(j));
    domain = string(factor_mapping.Domain(j));
    column_names(j) = factor_name + " (" + domain + ")";
end

question_table = array2table( ...
    question_matrix, ...
    'VariableNames',cellstr(column_names));

question_table.Properties.RowNames = ...
    cellstr("Item " + string((1:n_items)'));

end

function [average_table,interpretation_table] = ...
    create_layer2_interpretation_tables( ...
        B2_aligned,factor_mapping,original_A2_indices,zero_tol)

K1 = size(B2_aligned,1);
K2 = size(B2_aligned,2)-1;
n_runs = size(B2_aligned,3);

if numel(original_A2_indices) ~= K2
    error('Original Layer-2 index mapping does not match aligned B2.');
end

average_B2 = mean(B2_aligned,3);

factor_names = string(factor_mapping.DisplayFactor);
domains = string(factor_mapping.Domain);

% TRUE original Layer-2 names
layer2_names = "A2_" + string(original_A2_indices(:)');

average_table = array2table( ...
    average_B2(:,2:end), ...
    'VariableNames',cellstr(layer2_names));

average_table = addvars( ...
    average_table, ...
    factor_names, ...
    domains, ...
    average_B2(:,1), ...
    'Before',1, ...
    'NewVariableNames',{ ...
    'Factor', ...
    'Domain', ...
    'MeanIntercept'});

n_rows = K1*K2;

DisplayA2 = strings(n_rows,1);
OriginalA2 = strings(n_rows,1);
DisplayFactor = strings(n_rows,1);
OriginalFactor = strings(n_rows,1);
Domain = strings(n_rows,1);

MeanIntercept = zeros(n_rows,1);
MeanLoading = zeros(n_rows,1);
SDLoading = zeros(n_rows,1);

SelectionPercent = zeros(n_rows,1);
PositivePercent = zeros(n_rows,1);
NegativePercent = zeros(n_rows,1);
SignConsistencyPercent = zeros(n_rows,1);

MeanProbabilityChange = zeros(n_rows,1);
SDProbabilityChange = zeros(n_rows,1);

logistic = @(x) 1./(1+exp(-x));

row = 0;

for h = 1:K2

    for k = 1:K1

        row = row + 1;

        intercept_runs = reshape(B2_aligned(k,1,:),n_runs,1);
        loading_runs = reshape(B2_aligned(k,h+1,:),n_runs,1);

        nonzero = abs(loading_runs) > zero_tol;

        probability_change_runs = ...
            logistic(intercept_runs+loading_runs) - ...
            logistic(intercept_runs);

        DisplayA2(row) = layer2_names(h);
        OriginalA2(row) = layer2_names(h);

        DisplayFactor(row) = factor_names(k);
        OriginalFactor(row) = factor_names(k);
        Domain(row) = domains(k);

        MeanIntercept(row) = mean(intercept_runs);
        MeanLoading(row) = mean(loading_runs);
        SDLoading(row) = std(loading_runs);

        SelectionPercent(row) = 100*mean(nonzero);
        PositivePercent(row) = 100*mean(loading_runs > zero_tol);
        NegativePercent(row) = 100*mean(loading_runs < -zero_tol);

        if any(nonzero)

            positive_given_selected = mean(loading_runs(nonzero)>0);
            negative_given_selected = mean(loading_runs(nonzero)<0);

            SignConsistencyPercent(row) = ...
                100*max(positive_given_selected,negative_given_selected);
        end

        MeanProbabilityChange(row) = mean(probability_change_runs);
        SDProbabilityChange(row) = std(probability_change_runs);
    end
end

interpretation_table = table( ...
    DisplayA2, ...
    OriginalA2, ...
    DisplayFactor, ...
    OriginalFactor, ...
    Domain, ...
    MeanIntercept, ...
    MeanLoading, ...
    SDLoading, ...
    SelectionPercent, ...
    PositivePercent, ...
    NegativePercent, ...
    SignConsistencyPercent, ...
    MeanProbabilityChange, ...
    SDProbabilityChange);

interpretation_table.AbsMeanLoading = ...
    abs(interpretation_table.MeanLoading);

interpretation_table = sortrows( ...
    interpretation_table, ...
    {'DisplayA2','AbsMeanLoading'}, ...
    {'ascend','descend'});

end

%% ============================================================
% IPIP-FFM item lookup
%% ============================================================

function lookup = IPIP_item_lookup()

lookup = [ ...
    "EXT1","I am the life of the party."; ...
    "EXT2","I don't talk a lot."; ...
    "EXT3","I feel comfortable around people."; ...
    "EXT4","I keep in the background."; ...
    "EXT5","I start conversations."; ...
    "EXT6","I have little to say."; ...
    "EXT7","I talk to a lot of different people at parties."; ...
    "EXT8","I don't like to draw attention to myself."; ...
    "EXT9","I don't mind being the center of attention."; ...
    "EXT10","I am quiet around strangers."; ...

    "EST1","I get stressed out easily."; ...
    "EST2","I am relaxed most of the time."; ...
    "EST3","I worry about things."; ...
    "EST4","I seldom feel blue."; ...
    "EST5","I am easily disturbed."; ...
    "EST6","I get upset easily."; ...
    "EST7","I change my mood a lot."; ...
    "EST8","I have frequent mood swings."; ...
    "EST9","I get irritated easily."; ...
    "EST10","I often feel blue."; ...

    "AGR1","I feel little concern for others."; ...
    "AGR2","I am interested in people."; ...
    "AGR3","I insult people."; ...
    "AGR4","I sympathize with others' feelings."; ...
    "AGR5","I am not interested in other people's problems."; ...
    "AGR6","I have a soft heart."; ...
    "AGR7","I am not really interested in others."; ...
    "AGR8","I take time out for others."; ...
    "AGR9","I feel others' emotions."; ...
    "AGR10","I make people feel at ease."; ...

    "CSN1","I am always prepared."; ...
    "CSN2","I leave my belongings around."; ...
    "CSN3","I pay attention to details."; ...
    "CSN4","I make a mess of things."; ...
    "CSN5","I get chores done right away."; ...
    "CSN6","I often forget to put things back in their proper place."; ...
    "CSN7","I like order."; ...
    "CSN8","I shirk my duties."; ...
    "CSN9","I follow a schedule."; ...
    "CSN10","I am exacting in my work."; ...

    "OPN1","I have a rich vocabulary."; ...
    "OPN2","I have difficulty understanding abstract ideas."; ...
    "OPN3","I have a vivid imagination."; ...
    "OPN4","I am not interested in abstract ideas."; ...
    "OPN5","I have excellent ideas."; ...
    "OPN6","I do not have a good imagination."; ...
    "OPN7","I am quick to understand things."; ...
    "OPN8","I use difficult words."; ...
    "OPN9","I spend time reflecting on things."; ...
    "OPN10","I am full of ideas." ...
    ];

end

%% ============================================================
% Create Layer-2 domain summary
%% ============================================================

function summary_table = ...
    create_layer2_domain_summary(factor_cluster_table)

if isempty(factor_cluster_table)
    summary_table = table;
    return
end

A2_factors = unique( ...
    string(factor_cluster_table.DisplayA2), ...
    'stable');

summary_table = table;

for h = 1:numel(A2_factors)

    current_A2 = A2_factors(h);

    positive_rows = ...
        string(factor_cluster_table.DisplayA2)==current_A2 & ...
        string(factor_cluster_table.Direction)=="Positive";

    negative_rows = ...
        string(factor_cluster_table.DisplayA2)==current_A2 & ...
        string(factor_cluster_table.Direction)=="Negative";

    [positive_domains,positive_factors,positive_percent] = ...
        summarize_direction( ...
            factor_cluster_table(positive_rows,:));

    [negative_domains,negative_factors,negative_percent] = ...
        summarize_direction( ...
            factor_cluster_table(negative_rows,:));

    n_rows = max( ...
        numel(positive_domains), ...
        numel(negative_domains));

    if n_rows==0
        continue
    end

    positive_domains(end+1:n_rows,1) = "";
    positive_factors(end+1:n_rows,1) = "";
    positive_percent(end+1:n_rows,1) = NaN;

    negative_domains(end+1:n_rows,1) = "";
    negative_factors(end+1:n_rows,1) = "";
    negative_percent(end+1:n_rows,1) = NaN;

    positive_labels = strings(n_rows,1);
    negative_labels = strings(n_rows,1);

    for j = 1:n_rows

        if strlength(positive_domains(j))>0

            positive_labels(j) = compose( ...
                "%s (%.1f%%)", ...
                positive_domains(j), ...
                positive_percent(j));
        end

        if strlength(negative_domains(j))>0

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

    summary_table = [summary_table;current_summary];
end

end


%% ============================================================
% Summarize one Layer-2 direction
%% ============================================================

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

    idx = string(direction_table.Domain)==current_domain;

    domain_table = direction_table(idx,:);

    percentages(d) = mean(domain_table.SelectionPercent);

    [~,order] = sort( ...
        abs(domain_table.MeanLoading), ...
        'descend');

    domain_table = domain_table(order,:);

    factors(d) = strjoin( ...
        string(domain_table.DisplayFactor), ...
        ", ");
end

[percentages,order] = sort(percentages,'descend');

domains = domains(order);
factors = factors(order);

end


%% ============================================================
% Print Layer-2 domain summary
%% ============================================================

function print_layer2_domain_summary(summary_table)

if isempty(summary_table)

    fprintf('\nNo stable Layer-2 relationships were detected.\n');
    return
end

fprintf('\n');

fprintf('%-15s %-22s %-34s %-22s %-34s\n', ...
    'ActiveColumn', ...
    'PositiveDomain', ...
    'PositiveFactorsAssociated', ...
    'NegativeDomain', ...
    'NegativeFactorsAssociated');

separator = repmat('-',1,135);
fprintf('%s\n',separator);

previous_A2 = "";

for i = 1:height(summary_table)

    current_A2 = string(summary_table.ActiveColumn(i));

    if i>1 && current_A2~=previous_A2
        fprintf('%s\n',separator);
    end

    if i==1 || current_A2~=previous_A2
        active_label = current_A2;
    else
        active_label = "";
    end

    fprintf('%-15s %-22s %-34s %-22s %-34s\n', ...
        char(active_label), ...
        char(summary_table.PositiveDomain(i)), ...
        char(summary_table.PositiveFactorsAssociated(i)), ...
        char(summary_table.NegativeDomain(i)), ...
        char(summary_table.NegativeFactorsAssociated(i)));

    previous_A2 = current_A2;
end

fprintf('%s\n',separator);

end