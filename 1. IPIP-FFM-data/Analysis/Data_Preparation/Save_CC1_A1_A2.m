%% ============================================================
% Save canonical C1 A1 and A2 representations from all C1 fits
%
% Output:
%   27 fits x 16000 participants = 432000 rows
%
% Columns:
%   FitID
%   RepetitionID
%   ParticipantID
%   A1_1, ..., A1_12
%   A2_1, ..., A2_K2
%
% A1 and A2 are matched and reordered to correspond exactly
% to the post-processed IPIPFMM_CC1 B1 and B2.
%% ============================================================

clear
clc

%% Paths

root_dir = '/Users/avinandanroy/Desktop/Research/DDE/DDE_Fitting_Pshycology2';

fmm_dir = fullfile(root_dir,'1. IPIP-FFM-data');
fit_results_dir = fullfile(fmm_dir,'Fitting','Results');
postprocessed_dir = fullfile(fmm_dir,'Analysis','PostProcessed_Data');
predictive_dir = fullfile(fmm_dir,'Analysis','Predictive_Data');

addpath(fullfile(root_dir,'Utilities'));

if ~exist(predictive_dir,'dir')
    mkdir(predictive_dir);
end

results_file = fullfile(fit_results_dir,'IPIP_DDE_100_results6_comp.mat');
postprocessed_file = fullfile(postprocessed_dir,'IPIPFMM_CC1.mat');

output_file = fullfile(predictive_dir,'IPIPFMM_C1_A1_A2_27Fits.csv');

%% Load results and post-processed C1 object

load(results_file,'results');
load(postprocessed_file,'IPIPFMM_CC1');

repetition_ids = IPIPFMM_CC1.repetition_ids;

factor_order_original = IPIPFMM_CC1.factor_order_original;
active_A2_original = IPIPFMM_CC1.active_A2_original;

n_fits = numel(repetition_ids);

K1_final = numel(factor_order_original);
K2_final = numel(active_A2_original);

fprintf('Number of C1 fits: %d\n',n_fits);
fprintf('Number of final A1 factors: %d\n',K1_final);
fprintf('Number of final A2 factors: %d\n',K2_final);

%% Recover exact matching used in C1 post-processing

group_results = results(repetition_ids);

matched_fit = average_matched_loadings(group_results,0.1);

B1_assignments = matched_fit.B1_assignments;
B2_assignments = matched_fit.B2_assignments;

if size(B1_assignments,1) ~= n_fits
    error('Number of B1 assignments does not match number of C1 fits.');
end

if size(B2_assignments,1) ~= n_fits
    error('Number of B2 assignments does not match number of C1 fits.');
end

K1_full = size(B1_assignments,2);
K2_full = size(B2_assignments,2);

%% Determine participant count

A1_first = results(repetition_ids(1)).A1_new;

if size(A1_first,2) == K1_full

    N = size(A1_first,1);

elseif size(A1_first,1) == K1_full

    N = size(A1_first,2);

else

    error('Cannot determine orientation of A1_new.');

end

fprintf('Participants per fit: %d\n',N);
fprintf('Total output rows: %d\n',n_fits*N);

%% Allocate output

FitID = repelem((1:n_fits)',N);
RepetitionID = repelem(repetition_ids(:),N);
ParticipantID = repmat((1:N)',n_fits,1);

A1_all = zeros(n_fits*N,K1_final);
A2_all = zeros(n_fits*N,K2_final);

%% Match and reorder A1 and A2 for every C1 fit

for f = 1:n_fits

    r = repetition_ids(f);

    A1_current = results(r).A1_new;
    A2_current = results(r).A2_new;

    if isempty(A1_current)
        error('A1_new is empty for repetition %d.',r);
    end

    if isempty(A2_current)
        error('A2_new is empty for repetition %d.',r);
    end

    %% Orient A1 as N x K1

    if size(A1_current,2) == K1_full

        A1_current = double(A1_current);

    elseif size(A1_current,1) == K1_full

        A1_current = double(A1_current');

    else

        error('Unexpected A1_new dimensions in repetition %d.',r);

    end

    %% Orient A2 as N x K2

    if size(A2_current,2) == K2_full

        A2_current = double(A2_current);

    elseif size(A2_current,1) == K2_full

        A2_current = double(A2_current');

    else

        error('Unexpected A2_new dimensions in repetition %d.',r);

    end

    if size(A1_current,1) ~= N
        error('A1 participant count differs in repetition %d.',r);
    end

    if size(A2_current,1) ~= N
        error('A2 participant count differs in repetition %d.',r);
    end

    %% Match A1 using the same B1 assignment

    A1_matched = A1_current(:,B1_assignments(f,:));

    %% Apply final canonical C1 A1 ordering

    A1_final = A1_matched(:,factor_order_original);

    %% Match A2 using the same B2 assignment

    A2_matched = A2_current(:,B2_assignments(f,:));

    %% Keep final active C1 A2 factors

    A2_final = A2_matched(:,active_A2_original);

    %% Store

    row_idx = (f-1)*N + (1:N);

    A1_all(row_idx,:) = A1_final;
    A2_all(row_idx,:) = A2_final;

end

%% Construct output table

A1_names = cellstr("A1_" + string(1:K1_final));
A2_names = cellstr("A2_" + string(1:K2_final));

A1_table = array2table(A1_all,'VariableNames',A1_names);
A2_table = array2table(A2_all,'VariableNames',A2_names);

output_table = table( ...
    FitID, ...
    RepetitionID, ...
    ParticipantID);

output_table = [output_table,A1_table,A2_table];

%% Sanity checks

assert(height(output_table) == n_fits*N);
assert(width(A1_table) == K1_final);
assert(width(A2_table) == K2_final);

if ~all(ismember(A1_all(:),[0 1]))
    warning('Some A1 values are not binary 0/1.');
end

if ~all(ismember(A2_all(:),[0 1]))
    warning('Some A2 values are not binary 0/1.');
end

fprintf('\nFinal output dimensions: %d rows x %d columns\n', ...
    height(output_table),width(output_table));

%% Save

writetable(output_table,output_file);

fprintf('\nSaved canonical C1 A1/A2 representations to:\n%s\n',output_file);