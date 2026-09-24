clear
clc
close all

%% ============================================================
% IPIP-100 Post Processing
%
% Final representation:
%
% Survey rows:
% EXT | EST | AGR | CSN | OPN
%
% Layer-1 factors:
% EXT | EST | AGR | CSN | OPN
%
% Within each domain, survey items are ordered from highest
% to lowest loading on the first Layer-1 factor of that domain.
%
% After post-processing, Layer-1 factors are relabeled
% canonically as A1_1,...,A1_K.
%% ============================================================

%% Paths

root_dir = '/Users/avinandanroy/Desktop/Research/DDE/DDE_Fitting_Pshycology2';

ipip_dir = fullfile(root_dir,'2. IPIP 100');
fit_results_dir = fullfile(ipip_dir,'Fitting','Results');
postprocessing_dir = fullfile(ipip_dir,'PostProcessing');
analysis_dir = fullfile(ipip_dir,'Analysis','PostProcessed_Data');

addpath(fullfile(root_dir,'Utilities'));
addpath(postprocessing_dir);

if ~exist(analysis_dir,'dir')
    mkdir(analysis_dir);
end

load(fullfile(fit_results_dir,'IPIP100_DDE_fit.mat'),'results');

%% ============================================================
% 1. Original survey-question domain membership
%
% Original rows correspond to q1,...,q100.
%% ============================================================

EXT_idx = [3 10 14 18 23 29 33 39 43 49 53 59 63 69 73 79 83 89 93 99];

EST_idx = [11 12 17 19 27 30 37 40 47 50 57 60 67 70 77 80 87 90 97 100];

AGR_idx = [2 6 9 13 22 26 32 36 42 46 52 56 62 66 72 76 82 86 92 96];

CSN_idx = [5 8 15 20 25 28 35 38 45 48 55 58 65 68 75 78 85 88 95 98];

OPN_idx = [1 4 7 16 21 24 31 34 41 44 51 54 61 64 71 74 81 84 91 94];

domain_row_order = [EXT_idx EST_idx AGR_idx CSN_idx OPN_idx];

assert(numel(domain_row_order) == 100);
assert(isequal(sort(domain_row_order),1:100));

domain_names = ["EXT","EST","AGR","CSN","OPN"];
block_size = 20;
block_sizes = repmat(block_size,1,5);

%% ============================================================
% 2. Active Layer-1 factors and desired ordering
%
% Original fitted factors:
%
% EXT : A1_1,  A1_2
% EST : A1_8,  A1_18
% AGR : A1_13
% CSN : A1_5,  A1_32
% OPN : A1_10
%% ============================================================

groups_100 = {[1 2],[8 18],[13],[5 32],[10]};

factor_order = [groups_100{:}];

K1 = numel(factor_order);

active_A1 = find(any(abs(results.B1_CSP(:,2:end)) > 1e-12,1));

fprintf('Active A1 factors in original fit:\n');
disp(active_A1)

fprintf('Desired A1 ordering:\n');
disp(factor_order)

assert(isequal(sort(active_A1),sort(factor_order)), ...
    'groups_100 does not contain exactly the active A1 factors.');

%% ============================================================
% 3. Identify active Layer-2 factor
%% ============================================================

active_A2 = find(any(abs(results.B2_CSP(factor_order,2:end)) > 1e-12,1));

K2 = numel(active_A2);

fprintf('Active A2 factor(s) in original fit:\n');
disp(active_A2)

assert(K2 == 1,'Expected exactly one active A2 factor.');

%% ============================================================
% 4. First organize B1 rows and columns by domain
%
% Rows:
% EXT | EST | AGR | CSN | OPN
%
% Columns:
%
% Intercept |
% A1_1 A1_2 |
% A1_8 A1_18 |
% A1_13 |
% A1_5 A1_32 |
% A1_10
%% ============================================================

B1_domain = results.B1_CSP(domain_row_order,[1 factor_order+1]);

%% ============================================================
% 5. Determine the anchor factor used to order each row block
%
% First canonical factor in each domain:
%
% EXT -> canonical A1_1 -> original A1_1
% EST -> canonical A1_3 -> original A1_8
% AGR -> canonical A1_5 -> original A1_13
% CSN -> canonical A1_6 -> original A1_5
% OPN -> canonical A1_8 -> original A1_10
%% ============================================================

n_factors_domain = cellfun(@numel,groups_100);

first_factor_position = [1,1+cumsum(n_factors_domain(1:end-1))];

% +1 because column 1 of B1 is the intercept
anchor_cols = first_factor_position+1;

anchor_factor_original = factor_order(first_factor_position);

fprintf('\nWithin-domain sorting factors:\n');

for d = 1:5
    fprintf('%s: canonical A1_%d = original A1_%d\n', ...
        domain_names(d),first_factor_position(d),anchor_factor_original(d));
end

%% ============================================================
% 6. Sort survey questions within each Big-Five domain
%
% Highest loading -> lowest loading on that domain's first A1.
%% ============================================================

within_domain_order = zeros(1,100);

for d = 1:5

    block = (d-1)*block_size+1:d*block_size;

    loading = B1_domain(block,anchor_cols(d));

    [~,local_order] = sort(loading,'descend');

    within_domain_order(block) = block(local_order);

end

%% Convert to original q1,...,q100 indices

final_row_order = domain_row_order(within_domain_order);

assert(numel(final_row_order) == 100);
assert(isequal(sort(final_row_order),1:100));

%% ============================================================
% 7. Final canonical loading matrices
%% ============================================================

% 100 questions x (intercept + 8 active A1 factors)

B1 = results.B1_CSP(final_row_order,[1 factor_order+1]);

% 8 active A1 factors x (intercept + 1 active A2 factor)

B2 = results.B2_CSP(factor_order,[1 active_A2+1]);

%% ============================================================
% 8. Respondent-level latent variables
%
% A1 receives exactly the same permutation as the B1 columns
% and B2 rows.
%% ============================================================

A1 = results.A1_new(:,factor_order);

A2 = results.A2_new(:,active_A2);

%% ============================================================
% 9. Item-level latent Gaussian representation
%
% Same item permutation as the rows of B1.
%% ============================================================

Z = results.Z(:,final_row_order);

%% ============================================================
% 10. Canonical Layer-1 factor map
%
% After this point we refer only to:
%
% A1_1,...,A1_8
%
% Original fitted factor numbers are retained only for provenance.
%% ============================================================

canonical_factor = (1:K1)';
original_factor = factor_order';

canonical_name = "A1_" + string(canonical_factor);
original_name = "A1_" + string(original_factor);

factor_domain = strings(K1,1);

pos = 1;

for d = 1:5

    n_d = numel(groups_100{d});

    factor_domain(pos:pos+n_d-1) = domain_names(d);

    pos = pos+n_d;

end

factor_map = table( ...
    canonical_factor, ...
    original_factor, ...
    canonical_name, ...
    original_name, ...
    factor_domain, ...
    'VariableNames', ...
    {'CanonicalFactor','OriginalFactor','CanonicalName','OriginalName','Domain'});

%% ============================================================
% 11. Load exact q1,...,q100 question text
%% ============================================================

run(fullfile(postprocessing_dir,'code_map.m'));

original_variable = "q" + string((1:100)');

[found,idx] = ismember(original_variable,var_mapping.Variable);

assert(all(found),'Some questions could not be found in code_map.m.');

question_text = var_mapping.Question(idx);

%% ============================================================
% 12. Original Big-Five domain membership
%% ============================================================

original_domain = strings(100,1);

original_domain(EXT_idx) = "EXT";
original_domain(EST_idx) = "EST";
original_domain(AGR_idx) = "AGR";
original_domain(CSN_idx) = "CSN";
original_domain(OPN_idx) = "OPN";

assert(all(strlength(original_domain) > 0));

%% ============================================================
% 13. Final survey-item mapping
%
% Row j of item_map corresponds exactly to row j of B1
% and column j of Z.
%% ============================================================

canonical_row = (1:100)';
original_row = final_row_order';

variable = original_variable(final_row_order);
question = question_text(final_row_order);
domain = original_domain(final_row_order);

domain_row = repmat((1:block_size)',5,1);

item_map = table( ...
    canonical_row, ...
    original_row, ...
    variable, ...
    question, ...
    domain, ...
    domain_row, ...
    'VariableNames', ...
    {'CanonicalRow','OriginalRow','Variable','Question','Domain','DomainRow'});

%% ============================================================
% 14. Record the factor used for row sorting
%% ============================================================

sort_factor = strings(100,1);
sort_loading = zeros(100,1);

for d = 1:5

    block = (d-1)*block_size+1:d*block_size;

    sort_factor(block) = "A1_" + string(first_factor_position(d));

    sort_loading(block) = B1(block,anchor_cols(d));

end

item_map.SortFactor = sort_factor;
item_map.SortLoading = sort_loading;

%% ============================================================
% 15. Construct final compact IPIP-100 object
%% ============================================================

IPIP100 = struct();

IPIP100.dataset_name = "IPIP-100";

IPIP100.B1 = B1;
IPIP100.B2 = B2;

IPIP100.A1 = A1;
IPIP100.A2 = A2;

IPIP100.Z = Z;

IPIP100.item_map = item_map;
IPIP100.factor_map = factor_map;

IPIP100.domain_names = domain_names;
IPIP100.block_size = block_size;
IPIP100.block_sizes = block_sizes;

IPIP100.K1 = K1;
IPIP100.K2 = K2;

% Provenance
IPIP100.row_order = final_row_order;
IPIP100.factor_order_original = factor_order;
IPIP100.active_A2_original = active_A2;

%% ============================================================
% 16. Sanity checks
%% ============================================================

fprintf('\n============================================================\n');
fprintf('FINAL IPIP-100 REPRESENTATION\n');
fprintf('============================================================\n');

fprintf('B1 : %s\n',mat2str(size(IPIP100.B1)));
fprintf('B2 : %s\n',mat2str(size(IPIP100.B2)));
fprintf('A1 : %s\n',mat2str(size(IPIP100.A1)));
fprintf('A2 : %s\n',mat2str(size(IPIP100.A2)));
fprintf('Z  : %s\n',mat2str(size(IPIP100.Z)));

assert(size(B1,1) == 100);
assert(size(B1,2) == K1+1);

assert(size(B2,1) == K1);
assert(size(B2,2) == K2+1);

assert(size(A1,2) == K1);
assert(size(A2,2) == K2);

assert(size(Z,2) == 100);

assert(height(item_map) == 100);
assert(height(factor_map) == K1);

%% Check final domain blocks

expected_domains = [ ...
    repmat("EXT",block_size,1); ...
    repmat("EST",block_size,1); ...
    repmat("AGR",block_size,1); ...
    repmat("CSN",block_size,1); ...
    repmat("OPN",block_size,1)];

assert(isequal(item_map.Domain,expected_domains));

%% Check that rows are descending within every domain

for d = 1:5

    block = (d-1)*block_size+1:d*block_size;

    loading = B1(block,anchor_cols(d));

    assert(all(diff(loading) <= 1e-12), ...
        'Within-domain row ordering failed for %s.',domain_names(d));

end

fprintf('\nAll sanity checks passed.\n');

%% ============================================================
% 17. Save post-processed object and mapping files
%% ============================================================

save_file = fullfile(analysis_dir,'IPIP100_PostProcessed.mat');
item_map_file = fullfile(analysis_dir,'IPIP100_Item_Map.csv');
factor_map_file = fullfile(analysis_dir,'IPIP100_Factor_Map.csv');

save(save_file,'IPIP100','-v7.3');

writetable(item_map,item_map_file);
writetable(factor_map,factor_map_file);

%% ============================================================
% 18. Display mappings
%% ============================================================

fprintf('\n============================================================\n');
fprintf('CANONICAL LAYER-1 FACTORS\n');
fprintf('============================================================\n');

disp(factor_map)

fprintf('\n============================================================\n');
fprintf('FINAL SURVEY-ITEM ORDER\n');
fprintf('============================================================\n');

disp(item_map)

%% ============================================================
% 19. Plot canonical loading matrices
%
% plot_loadings2 performs NO further row or factor permutation.
%% ============================================================

plot_loadings2(IPIP100,block_size);

%% ============================================================
% Finished
%% ============================================================

fprintf('\n============================================================\n');
fprintf('IPIP-100 POST-PROCESSING COMPLETED\n');
fprintf('============================================================\n');

fprintf('\nPost-processed fit:\n%s\n',save_file);
fprintf('\nItem mapping:\n%s\n',item_map_file);
fprintf('\nFactor mapping:\n%s\n',factor_map_file);