%% ============================================================
% LOOPR BFI-2 DDE Hyperparameter Tuning using pMSE
%% ============================================================

clear
clc
close all

%% ============================================================
% Add folders
%% ============================================================

script_dir = fileparts(mfilename('fullpath'));
loopr_dir = fileparts(fileparts(script_dir));
root_dir = fileparts(loopr_dir);

data_dir = fullfile(loopr_dir,'Data');
results_dir = fullfile(loopr_dir,'Fitting','Results');

addpath(fullfile(root_dir,'Scripts'));
addpath(fullfile(root_dir,'Algorithms'));
addpath(fullfile(root_dir,'Utilities'));

rng(20260826);

%% ============================================================
% Load LOOPR data
%% ============================================================

data_file = fullfile(data_dir,'LOOPRDataAgeGender.csv');

LOOPR = readtable(data_file,'VariableNamingRule','preserve');

%% ============================================================
% Extract BFI1-BFI60 only
%% ============================================================

var_names = arrayfun(@(j) sprintf('BFI%d',j),1:60,'UniformOutput',false);

missing_vars = setdiff(var_names,LOOPR.Properties.VariableNames);

if ~isempty(missing_vars)
    error('The following BFI variables are missing: %s',strjoin(missing_vars,', '));
end

X = double(LOOPR{:,var_names});

[N,J] = size(X);

fprintf('Number of participants: %d\n',N);
fprintf('Number of BFI items: %d\n',J);

%% ============================================================
% Validate BFI data
%% ============================================================

if any(isnan(X),'all')
    error('LOOPR contains missing BFI responses.');
end

if any(~ismember(X,1:5),'all')
    error('Unexpected BFI response values. Expected only 1-5.');
end

fprintf('BFI data are complete and contain only values 1-5.\n');

%% ============================================================
% Construct rank matrix R
%% ============================================================

R = NaN(N,J);

for j = 1:J

    yj = X(:,j);

    [~,~,rank_codes] = unique(yj,'sorted');

    R(:,j) = rank_codes;

end

%% ============================================================
% Number of rank levels
%% ============================================================

Rlevels = zeros(1,J);

for j = 1:J
    Rlevels(j) = max(R(:,j));
end

%% ============================================================
% Initialize DDE ONCE
%% ============================================================

K1_max = floor(J/3);
K2_max = floor(K1_max/3);
epsilon = 0.00001;

fprintf('\nInitializing DDE...\n');
fprintf('K1_max = %d\n',K1_max);
fprintf('K2_max = %d\n',K2_max);

rng(20260826);

Z_init = simulate_gaussian_mixture(X,J);

[prop_in,B1_in,B2_in,gamma_in,A1_in,A2_in] = ...
    Normal_init_v2(Z_init,K1_max,K1_max,K2_max,epsilon);

%% ============================================================
% Hyperparameter grid
%% ============================================================

taus = [0.1,0.2,0.4];
t1s = [0.02,0.04,0.06];
t2s = [0.02,0.04,0.06];
temps = [0.7,0.8,0.9];

nComb = numel(taus)*numel(t1s)*numel(t2s)*numel(temps);

fprintf('\nTotal hyperparameter combinations: %d\n',nComb);

%% ============================================================
% pMSE settings
%% ============================================================

nRep = 50;
Kfold = 5;
nTrees = 100;
minLeafSize = 5;

% Real = 1
% Synthetic = 0

y_all = [ones(N,1);zeros(N,1)];
y_cat = categorical(y_all);

c = mean(y_all);       % = 0.5

%% ============================================================
% Result storage
%
% score_matrix columns:
%
% 1 avg_pMSE
% 2 sd_pMSE
% 3 se_pMSE
% 4 min_pMSE
% 5 max_pMSE
% 6 num_act1
% 7 num_act2
%% ============================================================

score_matrix = NaN(nComb,7);

% tau, t1, t2, temp + 7 scores
tuning_matrix = NaN(nComb,11);

% Store all repetition-level pMSE values
pMSE_matrix = NaN(nComb,nRep);

%% ============================================================
% Best-model storage
%% ============================================================

best_avg_pMSE = Inf;

best_settings = struct('tau',[],'t1',[],'t2',[],'temp',[]);

best_num_act1 = [];
best_num_act2 = [];
best_it_num = [];

best_B1_CSP = [];
best_B2_CSP = [];
best_prop_CSP = [];
best_gamma_CSP = [];
best_Z = [];

best_model_file = fullfile(results_dir,'LOOPR_best_DDE_tuning.mat');
checkpoint_file = fullfile(results_dir,'LOOPR_DDE_tuning_checkpoint.mat');

%% ============================================================
% Restore checkpoint if available
%% ============================================================

if isfile(checkpoint_file)

    C = load(checkpoint_file);

    if isfield(C,'score_matrix') && isequal(size(C.score_matrix),size(score_matrix))
        score_matrix = C.score_matrix;
    end

    if isfield(C,'tuning_matrix') && isequal(size(C.tuning_matrix),size(tuning_matrix))
        tuning_matrix = C.tuning_matrix;
    end

    if isfield(C,'pMSE_matrix') && isequal(size(C.pMSE_matrix),size(pMSE_matrix))
        pMSE_matrix = C.pMSE_matrix;
    end

    fprintf('\nLoaded tuning checkpoint.\n');
    fprintf('Completed combinations: %d of %d\n',sum(~isnan(score_matrix(:,1))),nComb);

end

%% ============================================================
% Restore current best model
%% ============================================================

if isfile(best_model_file)

    B = load(best_model_file);

    if isfield(B,'best_avg_pMSE'), best_avg_pMSE = B.best_avg_pMSE; end
    if isfield(B,'best_settings'), best_settings = B.best_settings; end
    if isfield(B,'best_num_act1'), best_num_act1 = B.best_num_act1; end
    if isfield(B,'best_num_act2'), best_num_act2 = B.best_num_act2; end
    if isfield(B,'best_it_num'), best_it_num = B.best_it_num; end
    if isfield(B,'best_B1_CSP'), best_B1_CSP = B.best_B1_CSP; end
    if isfield(B,'best_B2_CSP'), best_B2_CSP = B.best_B2_CSP; end
    if isfield(B,'best_prop_CSP'), best_prop_CSP = B.best_prop_CSP; end
    if isfield(B,'best_gamma_CSP'), best_gamma_CSP = B.best_gamma_CSP; end
    if isfield(B,'best_Z'), best_Z = B.best_Z; end

    fprintf('Current best average pMSE: %.8f\n',best_avg_pMSE);

end

%% ============================================================
% Start parallel pool
%% ============================================================

pool = gcp('nocreate');

if isempty(pool)
    pool = parpool('Processes');
end

pool.IdleTimeout = Inf;

fprintf('Parallel pool has %d workers.\n',pool.NumWorkers);

%% ============================================================
% Hyperparameter tuning
%% ============================================================

res_idx = 1;

for t = 1:numel(taus)

    for l_1 = 1:numel(t1s)

        for l_2 = 1:numel(t2s)

            for tt = 1:numel(temps)

                tau = taus(t);
                t1 = t1s(l_1);
                t2 = t2s(l_2);
                temp = temps(tt);

                %% Skip already completed combinations

                if ~isnan(score_matrix(res_idx,1))

                    fprintf(['Skipping combination %d of %d: ', ...
                        'tau=%.3f, t1=%.3f, t2=%.3f, temp=%.3f\n'], ...
                        res_idx,nComb,tau,t1,t2,temp);

                    res_idx = res_idx+1;

                    continue

                end

                fprintf('\n============================================================\n');
                fprintf('Combination %d of %d\n',res_idx,nComb);
                fprintf('tau=%.3f, t1=%.3f, t2=%.3f, temp=%.3f\n',tau,t1,t2,temp);
                fprintf('============================================================\n');

                %% ========================================================
                % Fit DDE
                %% ========================================================

                rng(20260826+res_idx);

                fit_timer = tic;

                [num_act1,num_act2,theta,theta2,gamma_CSP,p,q, ...
                    A1_sample_long,A2_sample_long,Z,prop_CSP, ...
                    B1_CSP,B2_CSP,A1_new,A2_new,it_num] = ...
                    get_SAEM_RL_CSP( ...
                    X,Z_init,R,Rlevels, ...
                    prop_in,gamma_in,B1_in,B2_in, ...
                    A1_in,A2_in,1,50, ...
                    t1,t2,temp,tau);

                fit_runtime = toc(fit_timer);

                fprintf('\nDDE fit completed.\n');
                fprintf('Active Layer-1 factors: %d\n',num_act1);
                fprintf('Active Layer-2 factors: %d\n',num_act2);
                fprintf('Iterations: %d\n',it_num);
                fprintf('Runtime: %.2f minutes\n',fit_runtime/60);

                %% ========================================================
                % pMSE evaluation
                %% ========================================================

                pMSE = NaN(nRep,1);

                parfor s = 1:nRep

                    %% Generate synthetic data

                    rng(20270000 + 1000*res_idx + s,'twister');

                    [X_sim,~] = generate_X_Cop_pred( ...
                        N,prop_CSP,B1_CSP,B2_CSP,gamma_CSP,X);

                    X_sim = double(X_sim);

                    if any(isnan(X_sim),'all')
                        error('Synthetic data contain missing values.');
                    end

                    if any(~ismember(X_sim,1:5),'all')
                        error('Synthetic data contain values outside 1-5.');
                    end

                    %% Pool real and synthetic observations

                    X_all = [X;X_sim];

                    %% 5-fold cross-validation

                    rng(20280000 + 1000*res_idx + s,'twister');

                    cv = cvpartition(y_all,'KFold',Kfold);

                    p_hat = zeros(2*N,1);

                    for k = 1:Kfold

                        idxTrain = training(cv,k);
                        idxTest = test(cv,k);

                        %% Random forest classifier

                        rf = TreeBagger( ...
                            nTrees, ...
                            X_all(idxTrain,:), ...
                            y_cat(idxTrain), ...
                            'Method','classification', ...
                            'OOBPrediction','Off', ...
                            'MinLeafSize',minLeafSize);

                        %% Held-out probabilities

                        [~,scores] = predict(rf,X_all(idxTest,:));

                        class_names = string(rf.ClassNames);

                        idx_real = find(class_names=="1",1);

                        if isempty(idx_real)
                            error('TreeBagger did not return the real-data class.');
                        end

                        p_hat(idxTest) = scores(:,idx_real);

                    end

                    %% pMSE

                    pMSE(s) = mean((p_hat-c).^2);

                end

                %% ========================================================
                % Summarize pMSE
                %% ========================================================

                avg_pMSE = mean(pMSE);
                sd_pMSE = std(pMSE);
                se_pMSE = sd_pMSE/sqrt(nRep);
                min_pMSE = min(pMSE);
                max_pMSE = max(pMSE);

                fprintf('\npMSE summary:\n');
                fprintf('Mean pMSE: %.8f\n',avg_pMSE);
                fprintf('SD pMSE: %.8f\n',sd_pMSE);
                fprintf('SE pMSE: %.8f\n',se_pMSE);
                fprintf('Min pMSE: %.8f\n',min_pMSE);
                fprintf('Max pMSE: %.8f\n',max_pMSE);

                %% ========================================================
                % Store result
                %% ========================================================

                pMSE_matrix(res_idx,:) = pMSE';

                score_matrix(res_idx,:) = [ ...
                    avg_pMSE, ...
                    sd_pMSE, ...
                    se_pMSE, ...
                    min_pMSE, ...
                    max_pMSE, ...
                    num_act1, ...
                    num_act2];

                tuning_matrix(res_idx,:) = [ ...
                    tau, ...
                    t1, ...
                    t2, ...
                    temp, ...
                    score_matrix(res_idx,:)];

                %% ========================================================
                % Update best model
                %% ========================================================

                if avg_pMSE < best_avg_pMSE

                    best_avg_pMSE = avg_pMSE;

                    best_settings = struct( ...
                        'tau',tau, ...
                        't1',t1, ...
                        't2',t2, ...
                        'temp',temp);

                    best_num_act1 = num_act1;
                    best_num_act2 = num_act2;
                    best_it_num = it_num;

                    best_B1_CSP = B1_CSP;
                    best_B2_CSP = B2_CSP;
                    best_prop_CSP = prop_CSP;
                    best_gamma_CSP = gamma_CSP;
                    best_Z = Z;

                    save( ...
                        best_model_file, ...
                        'best_settings', ...
                        'best_avg_pMSE', ...
                        'best_num_act1', ...
                        'best_num_act2', ...
                        'best_it_num', ...
                        'best_B1_CSP', ...
                        'best_B2_CSP', ...
                        'best_prop_CSP', ...
                        'best_gamma_CSP', ...
                        'best_Z', ...
                        'var_names', ...
                        '-v7.3');

                    fprintf('\n*** New best model ***\n');
                    fprintf('Best average pMSE: %.8f\n',best_avg_pMSE);

                end

                fprintf('Current best average pMSE: %.8f\n',best_avg_pMSE);

                %% ========================================================
                % Save checkpoint
                %% ========================================================

                completed_combinations = sum(~isnan(score_matrix(:,1)));

                save( ...
                    checkpoint_file, ...
                    'score_matrix', ...
                    'tuning_matrix', ...
                    'pMSE_matrix', ...
                    'completed_combinations', ...
                    'taus', ...
                    't1s', ...
                    't2s', ...
                    'temps', ...
                    'nRep', ...
                    'Kfold', ...
                    '-v7.3');

                %% Clear current fitted model

                clear theta theta2 gamma_CSP p q
                clear A1_sample_long A2_sample_long Z
                clear prop_CSP B1_CSP B2_CSP A1_new A2_new pMSE

                res_idx = res_idx+1;

            end

        end

    end

end

%% ============================================================
% Final tuning table
%% ============================================================

column_names = [ ...
    "tau", ...
    "t1", ...
    "t2", ...
    "temp", ...
    "avg_pMSE", ...
    "sd_pMSE", ...
    "se_pMSE", ...
    "min_pMSE", ...
    "max_pMSE", ...
    "num_act1", ...
    "num_act2"];

tuning_table = array2table( ...
    tuning_matrix, ...
    'VariableNames',cellstr(column_names));

tuning_table = sortrows(tuning_table,'avg_pMSE','ascend');

fprintf('\n============================================================\n');
fprintf('FINAL LOOPR TUNING RESULTS\n');
fprintf('============================================================\n');

disp(tuning_table);

%% ============================================================
% Save final tuning results
%% ============================================================

results_mat = fullfile(results_dir,'LOOPR_DDE_pMSE_tuning_results.mat');
results_csv = fullfile(results_dir,'LOOPR_DDE_pMSE_tuning_results.csv');

save( ...
    results_mat, ...
    'score_matrix', ...
    'tuning_matrix', ...
    'pMSE_matrix', ...
    'tuning_table', ...
    'taus', ...
    't1s', ...
    't2s', ...
    'temps', ...
    'nRep', ...
    'Kfold', ...
    '-v7.3');

writetable(tuning_table,results_csv);

%% ============================================================
% Best result
%% ============================================================

fprintf('\nBest settings:\n');
disp(best_settings);

fprintf('Best average pMSE: %.8f\n',best_avg_pMSE);
fprintf('Active Layer-1 factors: %d\n',best_num_act1);
fprintf('Active Layer-2 factors: %d\n',best_num_act2);