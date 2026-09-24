%% ============================================================
% Post-process Joe's original IPIP-98 DDE fit
%
% SOURCE:
%   2. IPIP 100/Fitting/Results/B5_analysis.mat
%
% QUESTION MAP:
%   2. IPIP 100/PostProcessing/code_map.m
%
% OUTPUTS:
%   IPIP98_PostProcessed.mat
%   IPIP98_Item_Map.csv
%   IPIP98_Factor_Map.csv
%% ============================================================

%clear
%clc
%close all

%% ============================================================
% Paths
%% ============================================================

root_dir = '/Users/avinandanroy/Desktop/Research/DDE/DDE_Fitting_Pshycology2';

ipip_dir = fullfile(root_dir,'2. IPIP 100');
fit_results_dir = fullfile(ipip_dir,'Fitting','Results');
postprocessing_dir = fullfile(ipip_dir,'PostProcessing');
analysis_dir = fullfile(ipip_dir,'Analysis','PostProcessed_Data');

% source the B5_analysis.mat from Joe's Github
source_file = fullfile(fit_results_dir,'B5_analysis.mat');

save_file = fullfile(analysis_dir,'IPIP98_PostProcessed.mat');
item_map_file = fullfile(analysis_dir,'IPIP98_Item_Map.csv');
factor_map_file = fullfile(analysis_dir,'IPIP98_Factor_Map.csv');

%% ============================================================
% Load Joe's original analysis
%% ============================================================

S = load(source_file);

%% ============================================================
% Extract fitted quantities
%% ============================================================

B1_full = S.best_B1_CSP;
B2_full = S.best_B2_CSP;

A1_full = S.best_A1;
A2_full = S.best_A2;

X_full = S.X;

var_names = string(S.var_names(:));

%% ============================================================
% Basic dimension checks
%% ============================================================

assert(isequal(size(X_full),[1500 98]));
assert(isequal(size(B1_full),[98 33]));
assert(isequal(size(B2_full),[32 11]));
assert(size(A1_full,1)==1500 && size(A1_full,2)==33);
assert(size(A2_full,1)==1500 && size(A2_full,2)==11);
assert(numel(var_names)==98);

%% ============================================================
% Recover original BIGFIVE item numbers
%% ============================================================

q_number = str2double(erase(var_names,"BIGFIVE_"));

assert(all(~isnan(q_number)));

expected_q = setdiff(1:100,[44 51],'stable')';

assert(isequal(q_number,expected_q), ...
    'The 98 variables do not match BIGFIVE_1:100 excluding 44 and 51.');

fprintf('\nVerified: BIGFIVE_44 and BIGFIVE_51 are excluded.\n');

%% ============================================================
% Load question text from code_map.m
%% ============================================================

run(fullfile(postprocessing_dir,'code_map.m'));

question_codes = "q" + string(q_number);

[found,location] = ismember(question_codes,var_mapping.Variable);

assert(all(found), ...
    'Some IPIP-98 variables could not be found in code_map.m.');

question_source = var_mapping.Question(location);

assert(numel(question_source)==98);

fprintf('All 98 retained items matched to question text in code_map.m.\n');

%% ============================================================
% Big Five domain definitions
%% ============================================================

EXT_idx = [3 10 14 18 23 29 33 39 43 49 53 59 63 69 73 79 83 89 93 99];

EST_idx = [11 12 17 19 27 30 37 40 47 50 57 60 67 70 77 80 87 90 97 100];

AGR_idx = [2 6 9 13 22 26 32 36 42 46 52 56 62 66 72 76 82 86 92 96];

CSN_idx = [5 8 15 20 25 28 35 38 45 48 55 58 65 68 75 78 85 88 95 98];

OPN_idx = [1 4 7 16 21 24 31 34 41 44 51 54 61 64 71 74 81 84 91 94];

%% Remove the two excluded political-voting items

OPN_idx = setdiff(OPN_idx,[44 51],'stable');

domain_names = ["EXT","EST","AGR","CSN","OPN"];

domain_items = { ...
    EXT_idx, ...
    EST_idx, ...
    AGR_idx, ...
    CSN_idx, ...
    OPN_idx};

block_sizes = cellfun(@numel,domain_items);

assert(isequal(block_sizes,[20 20 20 20 18]));

%% ============================================================
% Source rows for each domain
%% ============================================================

domain_rows = cell(1,5);
item_domain = strings(98,1);

for d = 1:5

    domain_rows{d} = find(ismember(q_number,domain_items{d}))';

    item_domain(domain_rows{d}) = domain_names(d);

end

assert(all(strlength(item_domain)>0));

%% ============================================================
% Initial domain row ordering
%
% EXT | EST | AGR | CSN | OPN
%% ============================================================

domain_row_order = [domain_rows{:}];

assert(numel(domain_row_order)==98);
assert(isequal(sort(domain_row_order),1:98));

%% ============================================================
% Active Layer-1 factors
%% ============================================================

active_A1 = find(any(abs(B1_full(:,2:end))>1e-12,1));

assert(numel(active_A1)==9, ...
    'Expected 9 active Layer-1 factors.');

fprintf('\nActive Layer-1 factors:\n');
disp(active_A1)

%% ============================================================
% Active Layer-2 factors
%% ============================================================

active_A2 = find(any(abs(B2_full(:,2:end))>1e-12,1));

assert(numel(active_A2)==1, ...
    'Expected 1 active Layer-2 factor.');

fprintf('Active Layer-2 factor:\n');
disp(active_A2)

%% ============================================================
% Domain strength of each active Layer-1 factor
%
% L2 norm of B1 loadings within each Big Five domain
%% ============================================================

domain_strength = zeros(numel(active_A1),5);

for k = 1:numel(active_A1)

    old_factor = active_A1(k);

    for d = 1:5

        rows = domain_rows{d};

        domain_strength(k,d) = norm(B1_full(rows,old_factor+1),2);

    end

end

%% ============================================================
% Assign each active factor to strongest domain
%% ============================================================

[~,domain_assignment] = max(domain_strength,[],2);

groups_original = cell(1,5);

for d = 1:5

    idx = find(domain_assignment==d);

    if isempty(idx)
        groups_original{d} = [];
        continue
    end

    strength_d = domain_strength(idx,d);

    [~,ord] = sort(strength_d,'descend');

    groups_original{d} = active_A1(idx(ord));

end

%% ============================================================
% Canonical Layer-1 ordering
%
% EXT | EST | AGR | CSN | OPN
%% ============================================================

factor_order = [groups_original{:}];

assert(numel(factor_order)==9);
assert(isequal(sort(factor_order),sort(active_A1)));

fprintf('\n============================================================\n');
fprintf('LAYER-1 DOMAIN ORGANIZATION\n');
fprintf('============================================================\n');

for d = 1:5
    fprintf('%s : [%s]\n',domain_names(d),num2str(groups_original{d}));
end

%% ============================================================
% Factor-domain metadata
%% ============================================================

factor_domains = strings(9,1);
factor_domain_strength = zeros(9,1);

pos = 1;

for d = 1:5

    group = groups_original{d};

    for j = 1:numel(group)

        old_factor = group(j);

        old_pos = find(active_A1==old_factor);

        factor_domains(pos) = domain_names(d);
        factor_domain_strength(pos) = domain_strength(old_pos,d);

        pos = pos+1;

    end

end

%% ============================================================
% Factor map
%% ============================================================

CanonicalFactor = (1:9)';
OriginalFactor = factor_order';

CanonicalName = "A1_" + string(CanonicalFactor);
OriginalName = "A1_" + string(OriginalFactor);

FactorDomain = factor_domains;
DomainStrength = factor_domain_strength;

factor_map = table( ...
    CanonicalFactor, ...
    OriginalFactor, ...
    CanonicalName, ...
    OriginalName, ...
    FactorDomain, ...
    DomainStrength);

%% ============================================================
% Reorder B1 first by domain
%% ============================================================

B1_domain = B1_full( ...
    domain_row_order, ...
    [1 factor_order+1]);

%% ============================================================
% Sort items within each domain
%
% Use first canonical Layer-1 factor belonging to the domain.
%% ============================================================

within_domain_order = 1:98;

start_row = 1;

for d = 1:5

    block = start_row:(start_row+block_sizes(d)-1);

    factor_idx = find(factor_domains==domain_names(d),1,'first');

    if ~isempty(factor_idx)

        loading = B1_domain(block,factor_idx+1);

        [~,local_order] = sort(loading,'descend');

        within_domain_order(block) = block(local_order);

    end

    start_row = start_row+block_sizes(d);

end

%% ============================================================
% Final survey row order
%% ============================================================

row_order = domain_row_order(within_domain_order);

assert(numel(row_order)==98);
assert(isequal(sort(row_order),1:98));

%% ============================================================
% Final B1
%% ============================================================

B1 = B1_full( ...
    row_order, ...
    [1 factor_order+1]);

%% ============================================================
% Final B2
%% ============================================================

B2 = B2_full( ...
    factor_order, ...
    [1 active_A2+1]);

%% ============================================================
% Final A1 and A2
%% ============================================================

A1 = A1_full(:,factor_order+1);
A2 = A2_full(:,active_A2+1);

%% ============================================================
% Final X
%% ============================================================

X = X_full(:,row_order);

%% ============================================================
% Item map
%% ============================================================

CanonicalRow = (1:98)';
SourceRow = row_order';

OriginalRow = q_number(row_order);
OriginalRow = OriginalRow(:);

Variable = var_names(row_order);
Variable = Variable(:);

Question = question_source(row_order);
Question = Question(:);

Domain = item_domain(row_order);
Domain = Domain(:);

%% ============================================================
% Meaningful item labels such as EXT1, EXT2, ...
%
% Numbering is based on the original retained-item ordering
% within each domain, not the post-plot sorting order.
%% ============================================================

Item = strings(98,1);

for d = 1:5

    rows = find(item_domain==domain_names(d));

    original_q = q_number(rows);

    [~,ord] = sort(original_q,'ascend');

    for j = 1:numel(ord)

        source_row = rows(ord(j));

        final_row = find(row_order==source_row);

        Item(final_row) = domain_names(d) + string(j);

    end

end

item_map = table( ...
    CanonicalRow, ...
    SourceRow, ...
    OriginalRow, ...
    Variable, ...
    Item, ...
    Question, ...
    Domain);

%% ============================================================
% Layer-2 map
%% ============================================================

CanonicalA2 = 1;
OriginalA2 = active_A2;

CanonicalName = "A2_1";
OriginalName = "A2_" + string(active_A2);

A2_map = table( ...
    CanonicalA2, ...
    OriginalA2, ...
    CanonicalName, ...
    OriginalName);

%% ============================================================
% Construct IPIP98 object
%% ============================================================

IPIP98 = struct();

IPIP98.dataset_name = "IPIP-98";

IPIP98.X = X;

IPIP98.B1 = B1;
IPIP98.B2 = B2;

IPIP98.A1 = A1;
IPIP98.A2 = A2;

IPIP98.item_map = item_map;
IPIP98.factor_map = factor_map;
IPIP98.A2_map = A2_map;

IPIP98.domain_names = domain_names;
IPIP98.block_sizes = block_sizes;

IPIP98.K1 = 9;
IPIP98.K2 = 1;

IPIP98.groups_original = groups_original;

IPIP98.factor_order_original = factor_order;

IPIP98.active_A1_original = active_A1;
IPIP98.active_A2_original = active_A2;

IPIP98.factor_domains = factor_domains;
IPIP98.domain_counts = cellfun(@numel,groups_original);

IPIP98.row_order = row_order;

IPIP98.excluded_items = [44 51];

%% ============================================================
% Keep tuning information
%% ============================================================

if isfield(S,'best_settings')
    IPIP98.best_settings = S.best_settings;
end

if isfield(S,'best_avg_pMSE')
    IPIP98.best_avg_pMSE = S.best_avg_pMSE;
end

%% ============================================================
% Sanity checks
%% ============================================================

assert(isequal(size(IPIP98.X),[1500 98]));

assert(isequal(size(IPIP98.B1),[98 10]));
assert(isequal(size(IPIP98.B2),[9 2]));

assert(isequal(size(IPIP98.A1),[1500 9]));
assert(isequal(size(IPIP98.A2),[1500 1]));

assert(height(IPIP98.item_map)==98);
assert(height(IPIP98.factor_map)==9);

assert(sum(IPIP98.block_sizes)==98);

assert(ismember('Question',IPIP98.item_map.Properties.VariableNames));
assert(all(strlength(IPIP98.item_map.Question)>0));

%% ============================================================
% Summary
%% ============================================================

fprintf('\n============================================================\n');
fprintf('IPIP-98 POST-PROCESSING COMPLETE\n');
fprintf('============================================================\n');

fprintf('X  : %s\n',mat2str(size(IPIP98.X)));
fprintf('B1 : %s\n',mat2str(size(IPIP98.B1)));
fprintf('B2 : %s\n',mat2str(size(IPIP98.B2)));
fprintf('A1 : %s\n',mat2str(size(IPIP98.A1)));
fprintf('A2 : %s\n',mat2str(size(IPIP98.A2)));

fprintf('\nLayer-1 factors by domain:\n');

for d = 1:5
    fprintf('%s : %d\n',domain_names(d),IPIP98.domain_counts(d));
end

fprintf('\nFactor map:\n');
disp(IPIP98.factor_map)

fprintf('\nFirst rows of item map:\n');
disp(IPIP98.item_map(1:min(10,height(IPIP98.item_map)),:))

%% ============================================================
% Save final object and maps
%% ============================================================

save(save_file,'IPIP98','-v7.3');

writetable(IPIP98.item_map,item_map_file);
writetable(IPIP98.factor_map,factor_map_file);

fprintf('\nSaved:\n');
fprintf('  %s\n',save_file);
fprintf('  %s\n',item_map_file);
fprintf('  %s\n',factor_map_file);