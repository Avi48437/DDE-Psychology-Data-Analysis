function anchor_map_table(M)

% ANCHOR_MAP_TABLE
%
% Table representation of the hierarchical anchor map.

%% ============================================================
% Header
%% ============================================================

fprintf('\n============================================================\n');
fprintf('%s\n',upper(M.dataset_name));
fprintf('============================================================\n');

%% ============================================================
% Layer-1 anchor items
%% ============================================================

fprintf('\n============================================================\n');
fprintf('LAYER-1 ANCHOR ITEMS\n');
fprintf('============================================================\n');

factors_per_block = 3;

n_blocks = ceil(M.K1/factors_per_block);

for b = 1:n_blocks

    first_col = (b-1)*factors_per_block+1;
    last_col = min(b*factors_per_block,M.K1);

    fprintf('\n============================================================\n');
    fprintf('Layer-1 factors: block %d of %d\n',b,n_blocks);
    fprintf('============================================================\n');

    disp(M.anchor_question_table(:,first_col:last_col))

end

%% ============================================================
% No Layer-2 factors
%% ============================================================

if M.K2==0

    fprintf('\nNo active Layer-2 factors.\n');

    return

end

%% ============================================================
% Layer-2 factor grouping
%% ============================================================

fprintf('\n============================================================\n');
fprintf('LAYER-2 FACTOR GROUPING\n');
fprintf('============================================================\n');

fprintf('%-15s %-20s %-32s %-20s %-32s\n', ...
    'ActiveColumn', ...
    'PositiveDomain', ...
    'PositiveFactorsAssociated', ...
    'NegativeDomain', ...
    'NegativeFactorsAssociated');

separator = repmat('-',1,125);

fprintf('%s\n',separator);

for i = 1:height(M.layer2_summary)

    fprintf('%-15s %-20s %-32s %-20s %-32s\n', ...
        char(M.layer2_summary.ActiveColumn(i)), ...
        char(M.layer2_summary.PositiveDomain(i)), ...
        char(M.layer2_summary.PositiveFactorsAssociated(i)), ...
        char(M.layer2_summary.NegativeDomain(i)), ...
        char(M.layer2_summary.NegativeFactorsAssociated(i)));

end

fprintf('%s\n',separator);

end