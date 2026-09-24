%% ============================================================
% LOOPR DDE Post-Processing
%
% Domain order:
% EXT -> NEM -> AGR -> CSN -> OPN
%% ============================================================

clear
clc
close all

%% ============================================================
% Paths
%% ============================================================

script_dir = fileparts(mfilename('fullpath'));
loopr_dir = fileparts(script_dir);
root_dir = fileparts(loopr_dir);
analysis_dir = fullfile(loopr_dir,'Analysis','PostProcessed_Data');

addpath(fullfile(root_dir,'Scripts'));
addpath(fullfile(root_dir,'Algorithms'));
addpath(fullfile(root_dir,'Utilities'));

if ~exist(analysis_dir,'dir')
    mkdir(analysis_dir);
end

fit_file = fullfile(loopr_dir,'Fitting','Results','LOOPR_DDE_fit.mat');
item_text_file = fullfile(analysis_dir,'LOOPR_BFI2_Item_Text.csv');
post_file = fullfile(analysis_dir,'LOOPR_PostProcessed.mat');

%% ============================================================
% Load fitted DDE
%% ============================================================

F = load(fit_file);

B1 = F.B1_CSP;
B2 = F.B2_CSP;

A1 = F.A1_new;
A2 = F.A2_new;

X = F.X;
Z = F.Z;

num_act1 = F.num_act1;
num_act2 = F.num_act2;

[N,J] = size(X);

fprintf('\nLoaded LOOPR DDE fit.\n');
fprintf('N  = %d\n',N);
fprintf('J  = %d\n',J);
fprintf('K1 = %d\n',num_act1);
fprintf('K2 = %d\n',num_act2);

%% ============================================================
% Load BFI-2 item text
%% ============================================================

if ~isfile(item_text_file)
    error('BFI-2 item-text file not found:\n%s',item_text_file);
end

BFI_Item_Text = readtable(item_text_file,'TextType','string');

required_columns = {'OriginalBFIItem','Item','Question'};

if ~all(ismember(required_columns,BFI_Item_Text.Properties.VariableNames))
    error('LOOPR_BFI2_Item_Text.csv must contain OriginalBFIItem, Item, and Question.');
end

if height(BFI_Item_Text) ~= 60
    error('Expected 60 BFI-2 items in LOOPR_BFI2_Item_Text.csv.');
end

if ~isequal(sort(BFI_Item_Text.OriginalBFIItem(:)),(1:60)')
    error('OriginalBFIItem must contain each item number from 1 to 60 exactly once.');
end

fprintf('\nLoaded BFI-2 question text for %d items.\n',height(BFI_Item_Text));

%% ============================================================
% Group BFI-2 items by domain
%
% Desired order:
%
% EXT : 1,6,11,...,56
% NEM : 4,9,14,...,59
% AGR : 2,7,12,...,57
% CSN : 3,8,13,...,58
% OPN : 5,10,15,...,60
%% ============================================================

survey_order = [ ...
    1:5:60, ...
    4:5:60, ...
    2:5:60, ...
    3:5:60, ...
    5:5:60];

domain_names = ["EXT","NEM","AGR","CSN","OPN"];

B1 = B1(survey_order,:);

%% ============================================================
% Separate intercepts and loading matrices
%% ============================================================

B1_intercept = B1(:,1);
B1_loadings = B1(:,2:end);

B2_intercept = B2(:,1);
B2_loadings = B2(:,2:end);

if size(B1_loadings,2) ~= size(B2,1)
    error('Number of B1 factors does not match number of B2 rows.');
end

%% ============================================================
% Identify active Layer-1 factors
%% ============================================================

B1_strength = vecnorm(B1_loadings,2,1);
[~,strength_order] = sort(B1_strength,'descend');

B1_active = strength_order(1:num_act1);

%% ============================================================
% Compute domain strength of every Layer-1 factor
%% ============================================================

domain_strength = zeros(5,size(B1_loadings,2));

for d = 1:5
    rows = (d-1)*12 + (1:12);
    domain_strength(d,:) = vecnorm(B1_loadings(rows,:),2,1);
end

%% ============================================================
% Assign each Layer-1 factor to its strongest domain
%% ============================================================

[~,anchor_domain] = max(domain_strength,[],1);

%% ============================================================
% Order active Layer-1 factors
%
% EXT -> NEM -> AGR -> CSN -> OPN
%
% Within each domain:
% strongest factor first
%% ============================================================

B1_column_order = [];
domain_counts = zeros(1,5);

for d = 1:5

    idx = B1_active(anchor_domain(B1_active)==d);

    if ~isempty(idx)
        [~,ord] = sort(domain_strength(d,idx),'descend');
        idx = idx(ord);
    end

    domain_counts(d) = numel(idx);
    B1_column_order = [B1_column_order,idx];

end

if length(B1_column_order) ~= num_act1
    error('Not all active Layer-1 factors were assigned to a domain.');
end

fprintf('\nLayer-1 factors per domain:\n');

for d = 1:5
    fprintf('%s: %d\n',domain_names(d),domain_counts(d));
end

%% ============================================================
% Construct reordered B1
%% ============================================================

B1_final = [B1_intercept,B1_loadings(:,B1_column_order)];

%% ============================================================
% Apply identical Layer-1 factor order to B2 rows
%% ============================================================

B2_intercept = B2_intercept(B1_column_order,:);
B2_loadings = B2_loadings(B1_column_order,:);

%% ============================================================
% Reorder survey items within each domain
%
% Use the first Layer-1 factor assigned to each domain as
% reference and sort its loadings from highest to lowest.
%% ============================================================

final_survey_order = survey_order;

for d = 1:5

    rows = (d-1)*12 + (1:12);

    if domain_counts(d)==0
        continue
    end

    reference_column = 2 + sum(domain_counts(1:d-1));
    reference_values = B1_final(rows,reference_column);

    [~,row_perm] = sort(reference_values,'descend');

    current_block = B1_final(rows,:);
    B1_final(rows,:) = current_block(row_perm,:);

    current_order = final_survey_order(rows);
    final_survey_order(rows) = current_order(row_perm);

end

%% ============================================================
% Verify within-domain row ordering
%% ============================================================

for d = 1:5

    if domain_counts(d)==0
        continue
    end

    rows = (d-1)*12 + (1:12);
    reference_column = 2 + sum(domain_counts(1:d-1));
    vals = B1_final(rows,reference_column);

    if any(diff(vals)>1e-12)
        error('Item ordering failed for domain %s.',domain_names(d));
    end

end

fprintf('\nDomain-wise item ordering verified.\n');

%% ============================================================
% Identify active Layer-2 factors
%% ============================================================

B2_strength = vecnorm(B2_loadings,2,1);
[~,B2_strength_order] = sort(B2_strength,'descend');

B2_column_order = B2_strength_order(1:num_act2);

%% ============================================================
% Construct reordered B2
%% ============================================================

B2_final = [B2_intercept,B2_loadings(:,B2_column_order)];

%% ============================================================
% Apply identical factor permutations to respondent coordinates
%% ============================================================

A1_final = A1(:,B1_column_order);
A2_final = A2(:,B2_column_order);

%% ============================================================
% Apply survey-item permutation to X and Z
%% ============================================================

X_final = X(:,final_survey_order);
Z_final = Z(:,final_survey_order);

%% ============================================================
% Item map
%% ============================================================

CanonicalRow = (1:J)';
OriginalBFIItem = final_survey_order';
Item = "BFI" + string(OriginalBFIItem);
Domain = repelem(domain_names,12)';

[found_item,item_lookup_idx] = ismember(OriginalBFIItem,BFI_Item_Text.OriginalBFIItem);

if ~all(found_item)
    error('Could not match all LOOPR items to LOOPR_BFI2_Item_Text.csv.');
end

Question = string(BFI_Item_Text.Question(item_lookup_idx));

if ismember('OriginalVariable',BFI_Item_Text.Properties.VariableNames)
    OriginalVariable = string(BFI_Item_Text.OriginalVariable(item_lookup_idx));
    item_map = table(CanonicalRow,OriginalBFIItem,Item,OriginalVariable,Question,Domain);
else
    item_map = table(CanonicalRow,OriginalBFIItem,Item,Question,Domain);
end

%% ============================================================
% Verify question-text alignment
%% ============================================================

if any(ismissing(Question) | strlength(strtrim(Question))==0)
    error('One or more BFI-2 items have missing question text.');
end

fprintf('\nBFI-2 question-text alignment verified.\n');

%% ============================================================
% Layer-1 factor map
%% ============================================================

CanonicalFactor = (1:num_act1)';
OriginalFactor = B1_column_order';
Factor = "A1_" + string(CanonicalFactor);
FactorDomain = domain_names(anchor_domain(B1_column_order))';

factor_map = table(CanonicalFactor,OriginalFactor,Factor,FactorDomain);

%% ============================================================
% Layer-2 factor map
%% ============================================================

CanonicalA2 = (1:num_act2)';
OriginalA2 = B2_column_order';
A2Factor = "A2_" + string(CanonicalA2);

A2_map = table(CanonicalA2,OriginalA2,A2Factor);

%% ============================================================
% Construct post-processed LOOPR object
%% ============================================================

LOOPR = struct;

LOOPR.dataset_name = "LOOPR BFI-2";

LOOPR.X = X_final;
LOOPR.Z = Z_final;

LOOPR.B1 = B1_final;
LOOPR.B2 = B2_final;

LOOPR.A1 = A1_final;
LOOPR.A2 = A2_final;

LOOPR.item_map = item_map;
LOOPR.factor_map = factor_map;
LOOPR.A2_map = A2_map;

LOOPR.domain_names = domain_names;
LOOPR.block_size = 12;
LOOPR.block_sizes = repmat(12,1,5);

LOOPR.K1 = num_act1;
LOOPR.K2 = num_act2;

LOOPR.row_order = final_survey_order;
LOOPR.factor_order_original = B1_column_order;
LOOPR.active_A2_original = B2_column_order;
LOOPR.domain_counts = domain_counts;

if isfield(F,'best_settings')
    LOOPR.settings = F.best_settings;
end

if isfield(F,'prop_CSP')
    LOOPR.prop = F.prop_CSP;
end

if isfield(F,'gamma_CSP')
    LOOPR.gamma = F.gamma_CSP;
end

%% ============================================================
% Save mapping CSV files
%% ============================================================

item_map_file = fullfile(analysis_dir,'LOOPR_Item_Map.csv');
factor_map_file = fullfile(analysis_dir,'LOOPR_Factor_Map.csv');

writetable(LOOPR.item_map,item_map_file);
writetable(LOOPR.factor_map,factor_map_file);

%% ============================================================
% Save post-processed object
%% ============================================================

save(post_file,'LOOPR','-v7.3');

%% ============================================================
% Summary
%% ============================================================

fprintf('\n============================================================\n');
fprintf('LOOPR POST-PROCESSING COMPLETE\n');
fprintf('============================================================\n');
fprintf('N  = %d\n',N);
fprintf('J  = %d\n',J);
fprintf('K1 = %d\n',LOOPR.K1);
fprintf('K2 = %d\n',LOOPR.K2);
fprintf('Domain order: EXT -> NEM -> AGR -> CSN -> OPN\n');
fprintf('Domain counts: [%s]\n',num2str(domain_counts));

fprintf('\nSaved object:\n%s\n',post_file);
fprintf('\nItem map:\n%s\n',item_map_file);
fprintf('\nFactor map:\n%s\n',factor_map_file);
fprintf('============================================================\n');

%% ============================================================
% Plot post-processed loadings
%% ============================================================

plot_loadings2(LOOPR,12);