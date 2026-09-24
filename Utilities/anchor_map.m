function M = anchor_map(Fit,top_N,option)

% ANCHOR_MAP
%
% Extracts hierarchical anchor information from one already
% post-processed DDE fit.
%
% INPUTS
% -------
% Fit    : post-processed DDE object containing
%          Fit.B1
%          Fit.B2
%          Fit.item_map
%          Fit.factor_map
%
% top_N  : number of Layer-1 anchor items
%
% option : optional representation name
%
%          'table'   -> anchor_map_table.m
%          'network' -> anchor_map_network.m
%
%          Any future option automatically calls
%          anchor_map_<option>.m
%
% OUTPUT
% ------
% M      : common anchor-map structure

if nargin < 2 || isempty(top_N)
    top_N = 3;
end

%% ============================================================
% Validate
%% ============================================================

required_fields = {'B1','B2','item_map','factor_map'};

for j = 1:numel(required_fields)

    if ~isfield(Fit,required_fields{j})
        error('Fit is missing field %s.',required_fields{j});
    end

end

B1 = Fit.B1;
B2 = Fit.B2;

J = size(B1,1);
K1 = size(B1,2)-1;
K2 = size(B2,2)-1;

if height(Fit.item_map) ~= J
    error('item_map rows do not match B1 rows.');
end

if height(Fit.factor_map) ~= K1
    error('factor_map rows do not match B1 factors.');
end

if size(B2,1) ~= K1
    error('B2 rows do not match Layer-1 factors.');
end

%% ============================================================
% Dataset name
%% ============================================================

if isfield(Fit,'dataset_name')
    dataset_name = string(Fit.dataset_name);
else
    dataset_name = "DDE Fit";
end

%% ============================================================
% Layer-1 factor names
%% ============================================================

if ismember('CanonicalName',Fit.factor_map.Properties.VariableNames)

    factor_names = string(Fit.factor_map.CanonicalName);

elseif ismember('Factor',Fit.factor_map.Properties.VariableNames)

    factor_names = string(Fit.factor_map.Factor);

else

    factor_names = "A1_" + string((1:K1)');

end

factor_names = factor_names(:);

%% ============================================================
% Layer-1 domains
%% ============================================================

if ismember('Domain',Fit.factor_map.Properties.VariableNames)

    factor_domains = string(Fit.factor_map.Domain);

elseif ismember('FactorDomain',Fit.factor_map.Properties.VariableNames)

    factor_domains = string(Fit.factor_map.FactorDomain);

else

    factor_domains = repmat("",K1,1);

end

factor_domains = factor_domains(:);

%% ============================================================
% Item domains
%% ============================================================

if ~ismember('Domain',Fit.item_map.Properties.VariableNames)
    error('item_map must contain Domain.');
end

item_domains = string(Fit.item_map.Domain);
item_domains = item_domains(:);

%% ============================================================
% Original item identities
%% ============================================================

if ismember('OriginalRow',Fit.item_map.Properties.VariableNames)

    original_item = Fit.item_map.OriginalRow;

elseif ismember('OriginalBFIItem',Fit.item_map.Properties.VariableNames)

    original_item = Fit.item_map.OriginalBFIItem;

else

    error('item_map must contain OriginalRow or OriginalBFIItem.');
end

original_item = original_item(:);

%% ============================================================
% Construct within-domain item labels
%
% Examples:
%
% IPIP100:
% q18 -> EXT4
%
% LOOPR:
% BFI21 -> EXT5
%% ============================================================

item_labels = strings(J,1);

unique_domains = unique(item_domains,'stable');

for d = 1:numel(unique_domains)

    domain = unique_domains(d);

    rows = find(item_domains==domain);

    [~,ord] = sort(original_item(rows),'ascend');

    for j = 1:numel(ord)

        row = rows(ord(j));

        item_labels(row) = domain + string(j);

    end

end

%% ============================================================
% Question text
%% ============================================================

if ismember('Question',Fit.item_map.Properties.VariableNames)

    question_text = string(Fit.item_map.Question);

else

    question_text = repmat("",J,1);

end

question_text = question_text(:);

%% ============================================================
% Find Layer-1 anchor items
%
% For item j and factor k:
%
% tmp_jk =
%
% max{0,min_{l != k}(B1_jk-B1_jl)}
%
% Therefore an anchor item must load more strongly on factor k
% than on every competing Layer-1 factor.
%% ============================================================

top_items = strings(top_N,K1);
top_questions = strings(top_N,K1);

anchor_strength = zeros(top_N,K1);
anchor_indices = zeros(top_N,K1);

loading_cols = 2:(K1+1);

for k = 1:K1

    k_col = loading_cols(k);

    other_cols = setdiff(loading_cols,k_col);

    tmp = zeros(J,1);

    for j = 1:J

        if isempty(other_cols)

            tmp(j) = max(0,B1(j,k_col));

        else

            diff_load = B1(j,k_col)-B1(j,other_cols);

            tmp(j) = max(0,min(diff_load));

        end

    end

    [sorted_tmp,sort_idx] = sort(tmp,'descend');

    keep = find(sorted_tmp>0);

    n_keep = min(top_N,numel(keep));

    if n_keep>0

        idx = sort_idx(1:n_keep);

        top_items(1:n_keep,k) = item_labels(idx);
        top_questions(1:n_keep,k) = question_text(idx);

        anchor_strength(1:n_keep,k) = sorted_tmp(1:n_keep);
        anchor_indices(1:n_keep,k) = idx;

    end

end

%% ============================================================
% Display factor names
%% ============================================================

display_names = factor_names;

for k = 1:K1

    if strlength(factor_domains(k))>0

        display_names(k) = ...
            factor_names(k) + "_" + factor_domains(k);

    end

end

%% ============================================================
% Anchor item table
%% ============================================================

anchor_table = array2table( ...
    top_items, ...
    'VariableNames',cellstr(display_names), ...
    'RowNames',cellstr("Item " + string(1:top_N)));

%% ============================================================
% Anchor question table
%% ============================================================

question_display = strings(top_N,K1);

for k = 1:K1

    for j = 1:top_N

        if strlength(top_items(j,k))==0
            continue
        end

        if strlength(top_questions(j,k))>0

            question_display(j,k) = ...
                top_questions(j,k) + ...
                " [" + top_items(j,k) + "]";

        else

            question_display(j,k) = top_items(j,k);

        end

    end

end

anchor_question_table = array2table( ...
    question_display, ...
    'VariableNames',cellstr(display_names), ...
    'RowNames',cellstr("Item " + string(1:top_N)));

%% ============================================================
% Detailed Layer-2 table
%% ============================================================

B2_table = table;

for h = 1:K2

    A2_name = "A2_" + string(h);

    current = table( ...
        repmat(A2_name,K1,1), ...
        factor_names, ...
        factor_domains, ...
        B2(:,1), ...
        B2(:,h+1), ...
        'VariableNames', ...
        {'A2','A1','Domain','Intercept','Loading'});

    B2_table = [B2_table;current];

end

if ~isempty(B2_table)

    B2_table.AbsLoading = abs(B2_table.Loading);

    B2_table = sortrows( ...
        B2_table, ...
        {'A2','AbsLoading'}, ...
        {'ascend','descend'});

end

%% ============================================================
% Layer-2 factor grouping
%% ============================================================

layer2_summary = table;

for h = 1:K2

    A2_name = "A2_" + string(h);

    loading = B2(:,h+1);

    positive_idx = loading>0;
    negative_idx = loading<0;

    positive_domains = unique(factor_domains(positive_idx),'stable');
    negative_domains = unique(factor_domains(negative_idx),'stable');

    n_rows = max(numel(positive_domains),numel(negative_domains));

    if n_rows==0
        continue
    end

    PositiveDomain = strings(n_rows,1);
    PositiveFactorsAssociated = strings(n_rows,1);

    NegativeDomain = strings(n_rows,1);
    NegativeFactorsAssociated = strings(n_rows,1);

    %% Positive side

    for d = 1:numel(positive_domains)

        domain = positive_domains(d);

        idx = positive_idx & factor_domains==domain;

        factor_idx = find(idx);

        [~,ord] = sort(abs(loading(factor_idx)),'descend');

        factor_idx = factor_idx(ord);

        PositiveDomain(d) = domain;

        PositiveFactorsAssociated(d) = ...
            strjoin(factor_names(factor_idx),", ");

    end

    %% Negative side

    for d = 1:numel(negative_domains)

        domain = negative_domains(d);

        idx = negative_idx & factor_domains==domain;

        factor_idx = find(idx);

        [~,ord] = sort(abs(loading(factor_idx)),'descend');

        factor_idx = factor_idx(ord);

        NegativeDomain(d) = domain;

        NegativeFactorsAssociated(d) = ...
            strjoin(factor_names(factor_idx),", ");

    end

    ActiveColumn = strings(n_rows,1);
    ActiveColumn(1) = A2_name;

    current_summary = table( ...
        ActiveColumn, ...
        PositiveDomain, ...
        PositiveFactorsAssociated, ...
        NegativeDomain, ...
        NegativeFactorsAssociated);

    layer2_summary = [layer2_summary;current_summary];

end

%% ============================================================
% Construct common map structure
%% ============================================================

M = struct();

M.dataset_name = dataset_name;

M.K1 = K1;
M.K2 = K2;
M.top_N = top_N;

M.B1 = B1;
M.B2 = B2;

M.factor_names = factor_names;
M.factor_domains = factor_domains;
M.display_names = display_names;

M.item_labels = item_labels;
M.item_domains = item_domains;
M.question_text = question_text;

M.anchor_items = top_items;
M.anchor_questions = top_questions;
M.anchor_strength = anchor_strength;
M.anchor_indices = anchor_indices;

M.anchor_table = anchor_table;
M.anchor_question_table = anchor_question_table;

M.B2_table = B2_table;
M.layer2_summary = layer2_summary;

%% ============================================================
% Optional representation
%
% Example:
%
% anchor_map(Fit,3,'table')
%       -> anchor_map_table(M)
%
% anchor_map(Fit,3,'network')
%       -> anchor_map_network(M)
%
% Any future helper automatically works:
%
% anchor_map(Fit,3,'tree')
%       -> anchor_map_tree(M)
%% ============================================================

if nargin >= 3 && ~isempty(option)

    helper_name = "anchor_map_" + lower(string(option));

    if exist(char(helper_name),'file') ~= 2
        error('Representation helper %s.m was not found.',helper_name);
    end

    feval(char(helper_name),M);

end

end