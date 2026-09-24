%% DDE hyperparameter tuning using fixed Feldman-Kowal benchmarks
clear
clc
close all

%% Add folders
root_dir = pwd;
data_dir = fullfile(root_dir,'1. IPIP-FFM-data','Mice_Becnch_Result');

addpath(fullfile(root_dir,'Scripts'));
addpath(fullfile(root_dir,'Algorithms'));
addpath(fullfile(root_dir,'Utilities'));
addpath(data_dir);

rng(20260807);

%% Load the exact N = 16000 tuning sample
tuning_file = fullfile(data_dir,'X_tuning.mat');
S = load(tuning_file,'X_tuning');
X_tuning = S.X_tuning;

var_names = X_tuning.Properties.VariableNames;
X_obs = double(X_tuning{:,:});
X_obs(X_obs == 0) = NaN;

[N,J] = size(X_obs);

fprintf('Loaded X_tuning: %d rows x %d items\n',N,J);

%% Missingness summary
missing_mask = isnan(X_obs);
total_missing = sum(missing_mask,'all');
missing_percent = 100*total_missing/numel(X_obs);
rows_with_missing = sum(any(missing_mask,2));
complete_rows = sum(~any(missing_mask,2));

fprintf('\nMissingness summary:\n');
fprintf('Total missing responses: %d out of %d\n',total_missing,numel(X_obs));
fprintf('Overall missing percentage: %.4f%%\n',missing_percent);
fprintf('Respondents with at least one missing item: %d out of %d (%.2f%%)\n',rows_with_missing,N,100*rows_with_missing/N);
fprintf('Complete respondents: %d out of %d (%.2f%%)\n',complete_rows,N,100*complete_rows/N);
fprintf('Median missing items per respondent: %.1f\n',median(sum(missing_mask,2)));
fprintf('Maximum missing items for one respondent: %d\n',max(sum(missing_mask,2)));

%% Load the fixed benchmark datasets
benchmark_file = fullfile(data_dir,'X_Mice_Benchmark.csv');
benchmark_table = readtable(benchmark_file,'VariableNamingRule','preserve');

benchmark_names = string(benchmark_table.Properties.VariableNames);

% Remove optional benchmark ID column
benchmark_id_col = find(strcmpi(benchmark_names,'BenchmarkID'),1);

if ~isempty(benchmark_id_col)
    benchmark_table(:,benchmark_id_col) = [];
    benchmark_names = string(benchmark_table.Properties.VariableNames);
end

%% Match benchmark columns to X_tuning
tuning_names = string(var_names);
[found,location] = ismember(tuning_names,benchmark_names);

if any(~found)
    missing_names = tuning_names(~found);
    error('Benchmark file is missing the following X_tuning variables: %s',strjoin(missing_names,', '));
end

benchmark_table = benchmark_table(:,location);
X_Mice_Benchmark = double(benchmark_table{:,:});

clear benchmark_table

%% Validate benchmark dimensions
if mod(size(X_Mice_Benchmark,1),N) ~= 0
    error('Number of benchmark rows (%d) is not an integer multiple of N=%d.',size(X_Mice_Benchmark,1),N);
end

nBench = size(X_Mice_Benchmark,1)/N;

if any(isnan(X_Mice_Benchmark),'all')
    error('X_Mice_Benchmark contains missing values.');
end

if any(~ismember(X_Mice_Benchmark,1:5),'all')
    error('X_Mice_Benchmark contains values outside the IPIP support 1:5.');
end

% Store compactly
X_Mice_Benchmark = uint8(X_Mice_Benchmark);

fprintf('\nLoaded %d benchmark datasets.\n',nBench);
fprintf('Each benchmark has dimension %d x %d.\n',N,J);

%% Split the benchmark datasets
X_complete = cell(nBench,1);

for b = 1:nBench
    X_complete{b} = X_Mice_Benchmark((b-1)*N+(1:N),:);
end

clear X_Mice_Benchmark

%% Classifier labels
% Class 1 = benchmark
% Class 0 = synthetic
y_all = [ones(N,1);zeros(N,1)];
y_cat = categorical(y_all);
c = mean(y_all);

%% Construct missing-aware rank matrix
R = NaN(N,J);

for j = 1:J
    observed_idx = ~missing_mask(:,j);
    observed_values = X_obs(observed_idx,j);

    if isempty(observed_values)
        error('Item %s has no observed responses.',var_names{j});
    end

    [~,~,rank_codes] = unique(observed_values,'sorted');
    R(observed_idx,j) = rank_codes;
end

Rlevels = zeros(1,J);

for j = 1:J
    if all(isnan(R(:,j)))
        Rlevels(j) = 0;
    else
        Rlevels(j) = max(R(:,j));
    end
end

%% Median imputation used only for initialization
X_init = X_obs;
item_medians = NaN(1,J);

for j = 1:J
    observed_values = sort(X_obs(~missing_mask(:,j),j));

    if isempty(observed_values)
        error('Cannot calculate a median for item %s.',var_names{j});
    end

    item_medians(j) = observed_values(ceil(numel(observed_values)/2));
    X_init(missing_mask(:,j),j) = item_medians(j);
end

%% Initialize the DDE once
K1_max = floor(J/3);
K2_max = floor(K1_max/3);
epsilon = 0.00001;

fprintf('\nInitializing DDE with K1_max=%d and K2_max=%d...\n',K1_max,K2_max);

rng(20260807);
Z_init = simulate_gaussian_mixture(X_init,J);
[prop_in,B1_in,B2_in,gamma_in,A1_in,A2_in] = Normal_init_v2(Z_init,K1_max,K1_max,K2_max,epsilon);

%% Hyperparameter grid
taus = 0.1;
t1s = 0.02;
t2s = [0.02,0.04,0.06];
temps = [0.7,0.8,0.9];

nRep = 30;
Kfold = 5;

nComb = numel(taus)*numel(t1s)*numel(t2s)*numel(temps);

%% Result storage
% score_matrix columns:
% avg_pMSE, se_pMSE, max_pMSE, min_pMSE, num_act2
score_matrix = NaN(nComb,5);

% tuning_matrix columns:
% tau, t1, t2, temp, avg_pMSE, se_pMSE, max_pMSE, min_pMSE, num_act2
tuning_matrix = NaN(nComb,9);

%% Best-model storage
best_avg_pMSE = Inf;
best_se_pMSE = [];
best_max_pMSE = [];
best_min_pMSE = [];
best_num_act2 = [];
best_settings = struct('tau',[],'t1',[],'t2',[],'temp',[]);
best_it_num = [];

best_B1_CSP = [];
best_B2_CSP = [];
best_prop_CSP = [];
best_gamma_CSP = [];
best_Z = [];

best_model_file = fullfile(data_dir,'best_DDE_Mice_benchmark_tuning.mat');
checkpoint_file = fullfile(data_dir,'Mice_benchmark_tuning_checkpoint.mat');

%% Restore checkpoint if one exists
if isfile(checkpoint_file)
    C = load(checkpoint_file);

    if isfield(C,'score_matrix') && isequal(size(C.score_matrix),size(score_matrix))
        score_matrix = C.score_matrix;
    end

    if isfield(C,'tuning_matrix') && isequal(size(C.tuning_matrix),size(tuning_matrix))
        tuning_matrix = C.tuning_matrix;
    end

    fprintf('\nLoaded existing tuning checkpoint.\n');
    fprintf('Completed combinations found: %d of %d\n',sum(~isnan(score_matrix(:,1))),nComb);
end

%% Restore current best model if one exists
if isfile(best_model_file)
    B = load(best_model_file);

    if isfield(B,'best_avg_pMSE'), best_avg_pMSE = B.best_avg_pMSE; end
    if isfield(B,'best_se_pMSE'), best_se_pMSE = B.best_se_pMSE; end
    if isfield(B,'best_max_pMSE'), best_max_pMSE = B.best_max_pMSE; end
    if isfield(B,'best_min_pMSE'), best_min_pMSE = B.best_min_pMSE; end
    if isfield(B,'best_num_act2'), best_num_act2 = B.best_num_act2; end
    if isfield(B,'best_settings'), best_settings = B.best_settings; end
    if isfield(B,'best_it_num'), best_it_num = B.best_it_num; end
    if isfield(B,'best_B1_CSP'), best_B1_CSP = B.best_B1_CSP; end
    if isfield(B,'best_B2_CSP'), best_B2_CSP = B.best_B2_CSP; end
    if isfield(B,'best_prop_CSP'), best_prop_CSP = B.best_prop_CSP; end
    if isfield(B,'best_gamma_CSP'), best_gamma_CSP = B.best_gamma_CSP; end
    if isfield(B,'best_Z'), best_Z = B.best_Z; end

    fprintf('Loaded current best model with average pMSE %.8f.\n',best_avg_pMSE);
end

%% Hyperparameter tuning
res_idx = 1;

for t = 1:numel(taus)

    for l_1 = 1:numel(t1s)

        for l_2 = 1:numel(t2s)

            for tt = 1:numel(temps)

                %% Current hyperparameter combination
                tau = taus(t);
                t1 = t1s(l_1);
                t2 = t2s(l_2);
                temp = temps(tt);

                %% Skip combinations already completed
                if ~isnan(score_matrix(res_idx,1))
                    fprintf('Skipping combination %d of %d: tau=%.3f, t1=%.3f, t2=%.3f, temp=%.3f\n',res_idx,nComb,tau,t1,t2,temp);
                    res_idx = res_idx+1;
                    continue
                end

                fprintf('\n============================================================\n');
                fprintf('Combination %d of %d\n',res_idx,nComb);
                fprintf('tau=%.3f, t1=%.3f, t2=%.3f, temp=%.3f\n',tau,t1,t2,temp);

                %% Fit missing-aware DDE
                rng(20260807+res_idx);

                [num_act1,num_act2,theta,theta2,gamma_CSP,p,q,A1_sample_long,A2_sample_long,Z,prop_CSP,B1_CSP,B2_CSP,A1_new,A2_new,it_num] = get_SAEM_RL_CSP(X_init,Z_init,R,Rlevels,prop_in,gamma_in,B1_in,B2_in,A1_in,A2_in,1,50,t1,t2,temp,tau);

                fprintf('Active Layer-1 dimensions: %d\n',num_act1);
                fprintf('Active Layer-2 dimensions: %d\n',num_act2);
                fprintf('Iterations used: %d\n',it_num);

                %% Evaluate the benchmark datasets sequentially
                benchmark_mean_pMSE = NaN(nBench,1);

                for b = 1:nBench

                    fprintf('\nBenchmark %d of %d\n',b,nBench);

                    X_comp = double(X_complete{b});
                    pMSE_rep = NaN(nRep,1);

                    tBench = tic;

                    for s = 1:nRep

                        %% Generate synthetic dataset using margin adjustment
                        rng(20270000+1000*b+s,'twister');

                        [X_sim,~,~] = generate_X_Cop_pred_MA(N,prop_CSP,B1_CSP,B2_CSP,gamma_CSP,X_obs,Z);
                        X_sim = double(X_sim);

                        if any(isnan(X_sim),'all')
                            error('Synthetic data contain missing values for benchmark %d, repetition %d.',b,s);
                        end

                        if any(~ismember(X_sim,1:5),'all')
                            error('Synthetic data contain values outside 1:5 for benchmark %d, repetition %d.',b,s);
                        end

                        %% Benchmark versus synthetic data
                        X_all = [X_comp;X_sim];

                        %% Reproducible 5-fold partition
                        rng(20280000+1000*b+s,'twister');

                        cv = cvpartition(y_all,'KFold',Kfold);
                        p_hat = zeros(2*N,1);

                        for k = 1:Kfold

                            idxTrain = training(cv,k);
                            idxTest = test(cv,k);

                            rf = TreeBagger(100,X_all(idxTrain,:),y_cat(idxTrain), ...
                            'Method','classification','OOBPrediction','Off', ...
                            'MinLeafSize',5,'Options',statset('UseParallel',true));
                            [~,scores] = predict(rf,X_all(idxTest,:));

                            class_names = string(rf.ClassNames);
                            idx_real = find(class_names=="1",1);

                            if isempty(idx_real)
                                error('TreeBagger did not return the benchmark class.');
                            end

                            p_hat(idxTest) = scores(:,idx_real);
                        end

                        pMSE_rep(s) = mean((p_hat-c).^2);
                    end

                    %% Reduce repetitions to one benchmark mean
                    benchmark_mean_pMSE(b) = mean(pMSE_rep);

                    fprintf('Benchmark %d mean pMSE: %.8f\n',b,benchmark_mean_pMSE(b));
                    fprintf('Benchmark %d pMSE SD: %.8f\n',b,std(pMSE_rep));
                    fprintf('Benchmark %d wall time: %.1f min\n',b,toc(tBench)/60);
                end

                %% Summarize this hyperparameter combination
                avg_pMSE = mean(benchmark_mean_pMSE);
                se_pMSE = std(benchmark_mean_pMSE)/sqrt(nBench);
                max_pMSE = max(benchmark_mean_pMSE);
                min_pMSE = min(benchmark_mean_pMSE);

                %% Store score quantities
                score_matrix(res_idx,:) = [avg_pMSE,se_pMSE,max_pMSE,min_pMSE,num_act2];

                %% Map scores to hyperparameter settings
                tuning_matrix(res_idx,:) = [tau,t1,t2,temp,score_matrix(res_idx,:)];

                fprintf('\nCombination summary:\n');
                fprintf('Average pMSE: %.8f\n',avg_pMSE);
                fprintf('Benchmark-level SE: %.8f\n',se_pMSE);
                fprintf('Maximum benchmark pMSE: %.8f\n',max_pMSE);
                fprintf('Minimum benchmark pMSE: %.8f\n',min_pMSE);
                fprintf('Active Layer-2 dimensions: %d\n',num_act2);

                %% Save the best fitted DDE model
                if avg_pMSE < best_avg_pMSE

                    best_avg_pMSE = avg_pMSE;
                    best_se_pMSE = se_pMSE;
                    best_max_pMSE = max_pMSE;
                    best_min_pMSE = min_pMSE;
                    best_num_act2 = num_act2;
                    best_settings = struct('tau',tau,'t1',t1,'t2',t2,'temp',temp);
                    best_it_num = it_num;

                    best_B1_CSP = B1_CSP;
                    best_B2_CSP = B2_CSP;
                    best_prop_CSP = prop_CSP;
                    best_gamma_CSP = gamma_CSP;
                    best_Z = Z;

                    save(best_model_file,'best_settings','best_avg_pMSE','best_se_pMSE','best_max_pMSE','best_min_pMSE','best_num_act2','best_it_num','best_B1_CSP','best_B2_CSP','best_prop_CSP','best_gamma_CSP','best_Z','var_names','nBench','nRep','-v7.3');

                    fprintf('New best average pMSE: %.8f\n',best_avg_pMSE);
                end

                fprintf('Current best average pMSE: %.8f\n',best_avg_pMSE);

                %% Save checkpoint after every completed hyperparameter combination
                completed_combinations = sum(~isnan(score_matrix(:,1)));

                save(checkpoint_file,'score_matrix','tuning_matrix','completed_combinations','taus','t1s','t2s','temps','nBench','nRep');

                %% Move to next combination
                res_idx = res_idx+1;

                %% Clear current fitted DDE quantities
                clear theta theta2 gamma_CSP p q A1_sample_long A2_sample_long Z
                clear prop_CSP B1_CSP B2_CSP A1_new A2_new
                clear benchmark_mean_pMSE pMSE_rep X_comp X_sim X_all
            end
        end
    end
end

%% Create final tuning table
column_names = ["tau","t1","t2","temp","avg_pMSE","se_pMSE","max_pMSE","min_pMSE","num_act2"];

tuning_table = array2table(tuning_matrix,'VariableNames',cellstr(column_names));
tuning_table = sortrows(tuning_table,'avg_pMSE','ascend');

fprintf('\nFinal tuning results:\n');
disp(tuning_table);

%% Save final tuning results
results_mat = fullfile(data_dir,'Mice_benchmark_pMSE_tuning_results.mat');
results_csv = fullfile(data_dir,'Mice_benchmark_pMSE_tuning_results.csv');

save(results_mat,'score_matrix','tuning_matrix','tuning_table','taus','t1s','t2s','temps','nBench','nRep');
writetable(tuning_table,results_csv);

%% Final best result
fprintf('\nBest settings found:\n');
disp(best_settings);

fprintf('Best average pMSE: %.8f\n',best_avg_pMSE);
fprintf('Best benchmark-level SE: %.8f\n',best_se_pMSE);
fprintf('Best maximum benchmark pMSE: %.8f\n',best_max_pMSE);
fprintf('Best minimum benchmark pMSE: %.8f\n',best_min_pMSE);
fprintf('Best active Layer-2 dimensions: %d\n',best_num_act2);

fprintf('Best-model file:\n%s\n',best_model_file);
fprintf('Tuning results:\n%s\n',results_csv);