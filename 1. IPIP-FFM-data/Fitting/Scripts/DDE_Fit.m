%% ============================================================
% IPIP-FFM: 100 Repeated Complete-Case DDE Fits
%% ============================================================

%% Paths

root_dir = '/Users/avinandanroy/Desktop/Research/DDE/DDE_Fitting_Pshycology2';

fmm_dir = fullfile(root_dir,'1. IPIP-FFM-data');
data_dir = fullfile(fmm_dir,'Data');
results_dir = fullfile(fmm_dir,'Fitting','Results');

addpath(fullfile(root_dir,'Scripts'));
addpath(fullfile(root_dir,'Algorithms'));
addpath(fullfile(root_dir,'Utilities'));
addpath(fullfile(fmm_dir,'Helpers'));

if ~exist(results_dir,'dir')
    mkdir(results_dir);
end

results_file = fullfile(results_dir,'IPIP_DDE_100_results6_comp.mat');

%% Prepare complete-case IPIP-FFM data

load(fullfile(data_dir,'B5_Red.mat'),'B5_Red');

item_names = ["EXT"+string(1:10), ...
              "EST"+string(1:10), ...
              "AGR"+string(1:10), ...
              "CSN"+string(1:10), ...
              "OPN"+string(1:10)];

X_items = B5_Red{:,item_names};

complete_rows = all(ismember(X_items,1:5),2);
B5CM = B5_Red(complete_rows,:);

% Collapse countries with fewer than 1000 complete respondents

country = categorical(B5CM.country);
country_names = categories(country);
country_counts = countcats(country);

large_countries = country_names(country_counts >= 1000);

country_new = string(B5CM.country);
country_new(~ismember(country_new,string(large_countries))) = "OTHER";

B5CM.country = categorical(country_new);

fprintf('Complete-case respondents: %d\n',height(B5CM));

%% Settings

number_of_repetitions = 100;
sample_size = 16000;

K1_max = 16;
K2_max = 5;

tau = 0.10;
t1 = 0.02;
t2 = 0.02;
temp = 0.90;

epsilon = 0.00001;

%% Load existing results or initialize a new results file

if exist(results_file,'file')

    saved_data = load(results_file);

    results = saved_data.results;
    repetition_seeds = saved_data.repetition_seeds;
    completed = saved_data.completed;
    failed = saved_data.failed;

    fprintf('Existing results loaded.\n');
    fprintf('Completed repetitions: %d of %d\n',sum(completed),number_of_repetitions);

else

    rng('shuffle');

    repetition_seeds = randi(2^31 - 1,number_of_repetitions,1);

    completed = false(number_of_repetitions,1);
    failed = false(number_of_repetitions,1);

    results = repmat(struct( ...
        'seed', [], ...
        'runtime', [], ...
        'num_act1', [], ...
        'num_act2', [], ...
        'it_num', [], ...
        'theta', [], ...
        'theta2', [], ...
        'gamma_CSP', [], ...
        'p', [], ...
        'q', [], ...
        'A1_sample', [], ...
        'A2_sample', [], ...
        'Z', [], ...
        'prop_CSP', [], ...
        'B1_CSP', [], ...
        'B2_CSP', [], ...
        'A1_new', [], ...
        'A2_new', [], ...
        'B1_in', [], ...
        'B2_in', [], ...
        'gamma_in', [], ...
        'prop_in', [], ...
        'error_message', ''),number_of_repetitions,1);

    save(results_file, ...
        'results','repetition_seeds','completed','failed', ...
        'number_of_repetitions','sample_size', ...
        'K1_max','K2_max','tau','t1','t2', ...
        'temp','epsilon','-v7.3');

    fprintf('New results file created.\n');

end

%% Repeated fitting

for repetition = 1:number_of_repetitions

    if completed(repetition)
        fprintf('Repetition %d already completed. Skipping.\n',repetition);
        continue
    end

    fprintf('\n========================================\n');
    fprintf('Starting repetition %d of %d\n',repetition,number_of_repetitions);
    fprintf('Seed: %d\n',repetition_seeds(repetition));
    fprintf('========================================\n');

    rng(repetition_seeds(repetition),'twister');

    repetition_timer = tic;

    try

        %% Draw stratified complete-case sample

        Xmf = rndIPIPFMM(B5CM,sample_size);

        X = Xmf{:,:};

        [Nmf,J] = size(X);

        if any(~ismember(X(:),1:5))
            error('The sampled data contain invalid or missing item responses.');
        end

        %% Construct integer rank matrix

        Rmf = zeros(Nmf,J);

        for j = 1:J
            [~,~,integer_ranks] = unique(X(:,j),'sorted');
            Rmf(:,j) = integer_ranks;
        end

        Rlevelsmf = max(Rmf,[],1);

        %% Initialize DDE

        warning('off','all');

        Z_initmf = simulate_gaussian_mixture(X,J);

        [prop_inmf,B1_inmf,B2_inmf,gamma_inmf, ...
            A1_inmf,A2_inmf] = ...
            Normal_init_v2(Z_initmf,K1_max,K1_max,K2_max,epsilon);

        warning('on','all');

        %% Fit DDE copula

        [num_act1mf,num_act2mf,thetamf,theta2mf, ...
            gamma_CSPmf,pmf,qmf,A1_samplemf, ...
            A2_samplemf,Zmf,prop_CSPmf,B1_CSPmf, ...
            B2_CSPmf,A1_newmf,A2_newmf,it_nummf] = ...
            get_SAEM_RL_CSP( ...
                X,Z_initmf,Rmf,Rlevelsmf, ...
                prop_inmf,gamma_inmf,B1_inmf,B2_inmf, ...
                A1_inmf,A2_inmf,1,50, ...
                t1,t2,temp,tau);

        repetition_runtime = toc(repetition_timer);

        %% Store repetition results

        results(repetition).seed = repetition_seeds(repetition);
        results(repetition).runtime = repetition_runtime;

        results(repetition).num_act1 = num_act1mf;
        results(repetition).num_act2 = num_act2mf;
        results(repetition).it_num = it_nummf;

        results(repetition).theta = thetamf;
        results(repetition).theta2 = theta2mf;
        results(repetition).gamma_CSP = gamma_CSPmf;
        results(repetition).p = pmf;
        results(repetition).q = qmf;

        results(repetition).A1_sample = A1_samplemf;
        results(repetition).A2_sample = A2_samplemf;
        results(repetition).Z = Zmf;

        results(repetition).prop_CSP = prop_CSPmf;
        results(repetition).B1_CSP = B1_CSPmf;
        results(repetition).B2_CSP = B2_CSPmf;

        results(repetition).A1_new = A1_newmf;
        results(repetition).A2_new = A2_newmf;

        results(repetition).B1_in = B1_inmf;
        results(repetition).B2_in = B2_inmf;
        results(repetition).gamma_in = gamma_inmf;
        results(repetition).prop_in = prop_inmf;

        results(repetition).error_message = '';

        completed(repetition) = true;
        failed(repetition) = false;

        %% Save checkpoint after each repetition

        save(results_file, ...
            'results','repetition_seeds','completed','failed', ...
            'number_of_repetitions','sample_size', ...
            'K1_max','K2_max','tau','t1','t2', ...
            'temp','epsilon','-v7.3');

        fprintf('Repetition %d completed.\n',repetition);
        fprintf('Active layer-1 nodes: %d\n',num_act1mf);
        fprintf('Active layer-2 nodes: %d\n',num_act2mf);
        fprintf('Runtime: %.2f minutes\n',repetition_runtime/60);

    catch fit_error

        repetition_runtime = toc(repetition_timer);

        failed(repetition) = true;

        results(repetition).seed = repetition_seeds(repetition);
        results(repetition).runtime = repetition_runtime;
        results(repetition).error_message = fit_error.message;

        save(results_file, ...
            'results','repetition_seeds','completed','failed', ...
            'number_of_repetitions','sample_size', ...
            'K1_max','K2_max','tau','t1','t2', ...
            'temp','epsilon','-v7.3');

        fprintf(2,'Repetition %d failed: %s\n',repetition,fit_error.message);

    end
end

fprintf('\nCompleted repetitions: %d of %d\n',sum(completed),number_of_repetitions);