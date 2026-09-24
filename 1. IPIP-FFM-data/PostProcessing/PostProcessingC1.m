%% ============================================================
% IPIP-FFM C1 Post-Processing
%
% Final domain order:
% EXT -> EST -> AGR -> CSN -> OPN
%
% Final Layer-1 organization:
% EXT : A1_1,  A1_14
% EST : A1_2,  A1_12
% AGR : A1_4,  A1_15
% CSN : A1_6,  A1_3
% OPN : A1_11, A1_7, A1_5, A1_9
%% ============================================================

clear
clc
close all

%% ============================================================
% Paths
%% ============================================================

root_dir = '/Users/avinandanroy/Desktop/Research/DDE/DDE_Fitting_Pshycology2';

fmm_dir = fullfile(root_dir,'1. IPIP-FFM-data');
data_dir = fullfile(fmm_dir,'Data');
fit_results_dir = fullfile(fmm_dir,'Fitting','Results');
analysis_dir = fullfile(fmm_dir,'Analysis','PostProcessed_Data');
helper_dir = fullfile(fmm_dir,'Helpers');

addpath(fullfile(root_dir,'Scripts'));
addpath(fullfile(root_dir,'Algorithms'));
addpath(fullfile(root_dir,'Utilities'));
addpath(helper_dir);

if ~exist(analysis_dir,'dir')
    mkdir(analysis_dir);
end

result_file = fullfile(fit_results_dir,'IPIP_DDE_100_results6_comp.mat');

load(result_file,'results','completed');

zero_tol = 1e-10;

%% ============================================================
% 1. Completed repetitions
%% ============================================================

valid_iterations = find(completed);

if isempty(valid_iterations)
    error('No completed repetitions were found.');
end

fprintf('\nCompleted repetitions: %d\n',numel(valid_iterations));

%% ============================================================
% 2. Active dimensions for every repetition
%% ============================================================

n_valid = numel(valid_iterations);

K1_active = zeros(n_valid,1);
K2_active = zeros(n_valid,1);

for i = 1:n_valid

    r = valid_iterations(i);

    B1_r = results(r).B1_CSP;
    B2_r = results(r).B2_CSP;

    K1_active(i) = sum(any(abs(B1_r(:,2:end))>zero_tol,1));
    K2_active(i) = sum(any(abs(B2_r(:,2:end))>zero_tol,1));

end

%% ============================================================
% 3. Active-dimension combinations
%% ============================================================

dimension_pairs = [K1_active,K2_active];

[unique_pairs,~,group_id] = unique(dimension_pairs,'rows');

frequency = accumarray(group_id,1);
percent = 100*frequency/n_valid;

dimension_summary = table( ...
    unique_pairs(:,1), ...
    unique_pairs(:,2), ...
    frequency, ...
    percent, ...
    'VariableNames', ...
    {'K1_Active','K2_Active','Count','Percent'});

dimension_summary = sortrows( ...
    dimension_summary, ...
    {'Count','K1_Active','K2_Active'}, ...
    {'descend','ascend','ascend'});

dimension_summary.CombinationID = ...
    "C" + string((1:height(dimension_summary))');

dimension_summary = movevars( ...
    dimension_summary,'CombinationID','Before',1);

disp(dimension_summary)

%% ============================================================
% 4. Select C1
%% ============================================================

target_K1 = dimension_summary.K1_Active(1);
target_K2 = dimension_summary.K2_Active(1);

repetition_ids = valid_iterations( ...
    K1_active==target_K1 & ...
    K2_active==target_K2);

fprintf('\n============================================================\n');
fprintf('C1 CONFIGURATION\n');
fprintf('============================================================\n');
fprintf('K1   = %d\n',target_K1);
fprintf('K2   = %d\n',target_K2);
fprintf('Runs = %d\n',numel(repetition_ids));

%% ============================================================
% 5. Average matched C1 loading matrices
%
% average_matched_loadings already performs the required
% matching/alignment across C1 repetitions while retaining
% the original factor-column identities.
%% ============================================================

average_fit_comb = average_matched_loadings( ...
    results(repetition_ids),0.1);

average_B1_comb = average_fit_comb.average_B1;
average_B2_comb = average_fit_comb.average_B2;

fprintf('\nAverage C1 matrices:\n');
fprintf('B1: %s\n',mat2str(size(average_B1_comb)));
fprintf('B2: %s\n',mat2str(size(average_B2_comb)));

%% ============================================================
% 6. Final Layer-1 factor organization
%% ============================================================

groups_FMM = { ...
    [1 14], ...        % EXT
    [2 12], ...        % EST
    [4 15], ...        % AGR
    [6 3], ...         % CSN
    [11 7 5 9]};       % OPN

domain_names = ["EXT","EST","AGR","CSN","OPN"];

factor_order_original = [groups_FMM{:}];
domain_counts = cellfun(@numel,groups_FMM);

%% Verify requested Layer-1 factors

active_average_A1 = find( ...
    any(abs(average_B1_comb(:,2:end))>zero_tol,1));

missing_A1 = factor_order_original( ...
    ~ismember(factor_order_original,active_average_A1));

if ~isempty(missing_A1)
    error('Missing Layer-1 factors: %s', ...
        strjoin("A1_"+string(missing_A1),', '));
end

fprintf('\nFinal original Layer-1 order:\n');
disp(factor_order_original)

%% ============================================================
% 7. Construct canonical B1
%% ============================================================

B1_final = [ ...
    average_B1_comb(:,1), ...
    average_B1_comb(:,factor_order_original+1)];

K1_final = numel(factor_order_original);

%% ============================================================
% 8. Identify active Layer-2 factors
%% ============================================================

active_A2_original = find( ...
    any(abs(average_B2_comb(:,2:end))>zero_tol,1));

if isempty(active_A2_original)

    K2_final = 0;

    B2_final = average_B2_comb( ...
        factor_order_original,1);

else

    K2_final = numel(active_A2_original);

    B2_final = [ ...
        average_B2_comb(factor_order_original,1), ...
        average_B2_comb( ...
            factor_order_original, ...
            active_A2_original+1)];

end

fprintf('Active original Layer-2 factors: [%s]\n', ...
    num2str(active_A2_original));

fprintf('\nCanonical dimensions:\n');
fprintf('B1: %s\n',mat2str(size(B1_final)));
fprintf('B2: %s\n',mat2str(size(B2_final)));

%% ============================================================
% 9. Layer-1 factor metadata
%% ============================================================

factor_domains = strings(K1_final,1);

pos = 1;

for d = 1:5

    n_d = domain_counts(d);

    factor_domains(pos:pos+n_d-1) = domain_names(d);

    pos = pos+n_d;

end

CanonicalFactor = (1:K1_final)';
OriginalFactor = factor_order_original';

CanonicalName = "A1_" + string(CanonicalFactor);
OriginalName = "A1_" + string(OriginalFactor);

factor_map = table( ...
    CanonicalFactor, ...
    OriginalFactor, ...
    CanonicalName, ...
    OriginalName, ...
    factor_domains, ...
    'VariableNames', ...
    {'CanonicalFactor','OriginalFactor', ...
     'CanonicalName','OriginalName','Domain'});

%% ============================================================
% 10. Original IPIP-FFM item metadata
%% ============================================================

survey_vars = [ ...
    compose("EXT%d",1:10), ...
    compose("EST%d",1:10), ...
    compose("AGR%d",1:10), ...
    compose("CSN%d",1:10), ...
    compose("OPN%d",1:10)];

question_text = [ ...
    "I am the life of the party."
    "I don't talk a lot."
    "I feel comfortable around people."
    "I keep in the background."
    "I start conversations."
    "I have little to say."
    "I talk to a lot of different people at parties."
    "I don't like to draw attention to myself."
    "I don't mind being the center of attention."
    "I am quiet around strangers."

    "I get stressed out easily."
    "I am relaxed most of the time."
    "I worry about things."
    "I seldom feel blue."
    "I am easily disturbed."
    "I get upset easily."
    "I change my mood a lot."
    "I have frequent mood swings."
    "I get irritated easily."
    "I often feel blue."

    "I feel little concern for others."
    "I am interested in people."
    "I insult people."
    "I sympathize with others' feelings."
    "I am not interested in other people's problems."
    "I have a soft heart."
    "I am not really interested in others."
    "I take time out for others."
    "I feel others' emotions."
    "I make people feel at ease."

    "I am always prepared."
    "I leave my belongings around."
    "I pay attention to details."
    "I make a mess of things."
    "I get chores done right away."
    "I often forget to put things back in their proper place."
    "I like order."
    "I shirk my duties."
    "I follow a schedule."
    "I am exacting in my work."

    "I have a rich vocabulary."
    "I have difficulty understanding abstract ideas."
    "I have a vivid imagination."
    "I am not interested in abstract ideas."
    "I have excellent ideas."
    "I do not have a good imagination."
    "I am quick to understand things."
    "I use difficult words."
    "I spend time reflecting on things."
    "I am full of ideas."
    ];

item_domains = [ ...
    repmat("EXT",10,1); ...
    repmat("EST",10,1); ...
    repmat("AGR",10,1); ...
    repmat("CSN",10,1); ...
    repmat("OPN",10,1)];

item_map_original = table( ...
    (1:50)', ...
    survey_vars', ...
    question_text, ...
    item_domains, ...
    'VariableNames', ...
    {'OriginalRow','Variable','Question','Domain'});

%% ============================================================
% 11. Sort items within each domain
%
% First canonical Layer-1 factor in each domain is reference.
%% ============================================================

block_sizes = [10 10 10 10 10];

first_factor_position = ...
    [1,1+cumsum(domain_counts(1:end-1))];

reference_columns = first_factor_position+1;

final_row_order = 1:50;

for d = 1:5

    rows = (d-1)*10 + (1:10);

    reference_values = B1_final(rows,reference_columns(d));

    [~,row_perm] = sort(reference_values,'descend');

    current_block = B1_final(rows,:);
    B1_final(rows,:) = current_block(row_perm,:);

    current_order = final_row_order(rows);
    final_row_order(rows) = current_order(row_perm);

end

%% ============================================================
% 12. Verify within-domain item ordering
%% ============================================================

for d = 1:5

    rows = (d-1)*10 + (1:10);

    vals = B1_final(rows,reference_columns(d));

    if any(diff(vals)>1e-12)
        error('Within-domain ordering failed for %s.',domain_names(d));
    end

end

fprintf('\nWithin-domain item ordering verified.\n');

%% ============================================================
% 13. Final item map
%% ============================================================

item_map = item_map_original(final_row_order,:);

CanonicalRow = (1:50)';

item_map = addvars( ...
    item_map, ...
    CanonicalRow, ...
    'Before',1, ...
    'NewVariableNames','CanonicalRow');

DomainRow = repmat((1:10)',5,1);

item_map = addvars( ...
    item_map, ...
    DomainRow, ...
    'After','Domain', ...
    'NewVariableNames','DomainRow');

SortFactorCanonical = strings(50,1);
SortFactorOriginal = strings(50,1);
SortLoading = zeros(50,1);

for d = 1:5

    rows = (d-1)*10 + (1:10);

    SortFactorCanonical(rows) = ...
        "A1_" + string(first_factor_position(d));

    SortFactorOriginal(rows) = ...
        "A1_" + string( ...
        factor_order_original(first_factor_position(d)));

    SortLoading(rows) = ...
        B1_final(rows,reference_columns(d));

end

item_map.SortFactorCanonical = SortFactorCanonical;
item_map.SortFactorOriginal = SortFactorOriginal;
item_map.SortLoading = SortLoading;

%% ============================================================
% 14. Layer-2 factor map
%% ============================================================

CanonicalA2 = (1:K2_final)';
OriginalA2 = active_A2_original(:);

CanonicalA2Name = "A2_" + string(CanonicalA2);
OriginalA2Name = "A2_" + string(OriginalA2);

A2_map = table( ...
    CanonicalA2, ...
    OriginalA2, ...
    CanonicalA2Name, ...
    OriginalA2Name);

%% ============================================================
% 15. Construct canonical IPIP-FFM C1 object
%% ============================================================

IPIPFMM_CC1 = struct();

IPIPFMM_CC1.dataset_name = "IPIP-FFM C1";

IPIPFMM_CC1.B1 = B1_final;
IPIPFMM_CC1.B2 = B2_final;

IPIPFMM_CC1.item_map = item_map;
IPIPFMM_CC1.factor_map = factor_map;
IPIPFMM_CC1.A2_map = A2_map;

IPIPFMM_CC1.domain_names = domain_names;

IPIPFMM_CC1.block_size = 10;
IPIPFMM_CC1.block_sizes = block_sizes;

IPIPFMM_CC1.K1 = K1_final;
IPIPFMM_CC1.K2 = K2_final;

IPIPFMM_CC1.C1_K1 = target_K1;
IPIPFMM_CC1.C1_K2 = target_K2;

IPIPFMM_CC1.groups_original = groups_FMM;
IPIPFMM_CC1.factor_order_original = factor_order_original;
IPIPFMM_CC1.active_A2_original = active_A2_original;

IPIPFMM_CC1.factor_domains = factor_domains;
IPIPFMM_CC1.row_order = final_row_order;

IPIPFMM_CC1.repetition_ids = repetition_ids;
IPIPFMM_CC1.n_runs = numel(repetition_ids);

%% ============================================================
% 16. Sanity checks
%% ============================================================

assert(size(IPIPFMM_CC1.B1,1)==50);
assert(size(IPIPFMM_CC1.B1,2)==K1_final+1);

assert(size(IPIPFMM_CC1.B2,1)==K1_final);
assert(size(IPIPFMM_CC1.B2,2)==K2_final+1);

assert(height(IPIPFMM_CC1.item_map)==50);
assert(height(IPIPFMM_CC1.factor_map)==K1_final);

assert(isequal( ...
    IPIPFMM_CC1.factor_map.OriginalFactor, ...
    factor_order_original'));

expected_domains = [ ...
    repmat("EXT",10,1); ...
    repmat("EST",10,1); ...
    repmat("AGR",10,1); ...
    repmat("CSN",10,1); ...
    repmat("OPN",10,1)];

assert(isequal(IPIPFMM_CC1.item_map.Domain,expected_domains));

fprintf('\nAll IPIP-FFM C1 post-processing checks passed.\n');

%% ============================================================
% 17. Save mapping CSV files
%% ============================================================

item_map_file = fullfile( ...
    analysis_dir,'IPIPFMM_CC1_Item_Map.csv');

factor_map_file = fullfile( ...
    analysis_dir,'IPIPFMM_CC1_Factor_Map.csv');

writetable(IPIPFMM_CC1.item_map,item_map_file);
writetable(IPIPFMM_CC1.factor_map,factor_map_file);

%% ============================================================
% 18. Save canonical object
%% ============================================================

save_file = fullfile( ...
    analysis_dir,'IPIPFMM_CC1.mat');

save(save_file,'IPIPFMM_CC1','-v7.3');

%% ============================================================
% 19. Display final factor map
%% ============================================================

fprintf('\n============================================================\n');
fprintf('FINAL IPIP-FFM C1 FACTOR MAP\n');
fprintf('============================================================\n');

disp(IPIPFMM_CC1.factor_map)

%% ============================================================
% 20. Display final item map
%% ============================================================

fprintf('\n============================================================\n');
fprintf('FINAL IPIP-FFM C1 ITEM MAP\n');
fprintf('============================================================\n');

disp(IPIPFMM_CC1.item_map(:, ...
    {'CanonicalRow', ...
     'OriginalRow', ...
     'Variable', ...
     'Question', ...
     'Domain', ...
     'DomainRow', ...
     'SortFactorOriginal', ...
     'SortLoading'}));

%% ============================================================
% 21. Plot final post-processed loadings
%% ============================================================

plot_loadings2(IPIPFMM_CC1,10);

%% ============================================================
% Finished
%% ============================================================

fprintf('\n============================================================\n');
fprintf('IPIP-FFM C1 POST-PROCESSING COMPLETED\n');
fprintf('============================================================\n');

fprintf('\nSaved object:\n%s\n',save_file);
fprintf('\nItem map:\n%s\n',item_map_file);
fprintf('\nFactor map:\n%s\n',factor_map_file);
fprintf('============================================================\n');