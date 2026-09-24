function [pMSE_rep,mean_pMSE,sd_pMSE] = run_pmse_complete(X_complete,nRep,benchmark_id,prop_CSP,B1_CSP,B2_CSP,gamma_CSP,Kfold)

N = size(X_complete,1);

% Class 1 = benchmark, Class 0 = synthetic
y_all = [ones(N,1);zeros(N,1)];
y_cat = categorical(y_all);
c = mean(y_all);

pMSE_rep = NaN(nRep,1);

% Memory-control settings
batchSize = 5;
nWorkers = 2;

for batch_start = 1:batchSize:nRep
    batch_end = min(batch_start+batchSize-1,nRep);
    reps = batch_start:batch_end;
    n_batch = numel(reps);

    fprintf('  Repetitions %d-%d of %d\n',batch_start,batch_end,nRep);

    % Always start this batch with completely fresh workers
    delete(gcp('nocreate'));
    parpool('Processes',nWorkers);

    pMSE_batch = NaN(n_batch,1);

    parfor ii = 1:n_batch
        s = reps(ii);

        rng(20270000+1000*benchmark_id+s,'twister');

        [X_sim,~] = generate_X_Cop_pred(N,prop_CSP,B1_CSP,B2_CSP,gamma_CSP,X_complete);
        X_sim = double(X_sim);

        if any(isnan(X_sim),'all')
            error('Synthetic data contain missing values for benchmark %d, repetition %d.',benchmark_id,s);
        end

        if any(~ismember(X_sim,1:5),'all')
            error('Synthetic data contain values outside 1:5 for benchmark %d, repetition %d.',benchmark_id,s);
        end

        X_all = [X_complete;X_sim];

        rng(20280000+1000*benchmark_id+s,'twister');
        cv = cvpartition(y_all,'KFold',Kfold);
        p_hat = zeros(2*N,1);

        for k = 1:Kfold
            idxTrain = training(cv,k);
            idxTest = test(cv,k);

            rf = TreeBagger(100,X_all(idxTrain,:),y_cat(idxTrain),'Method','classification','OOBPrediction','Off','MinLeafSize',5);
            [~,scores] = predict(rf,X_all(idxTest,:));

            class_names = string(rf.ClassNames);
            idx_real = find(class_names=="1",1);

            if isempty(idx_real)
                error('TreeBagger did not return the benchmark class.');
            end

            p_hat(idxTest) = scores(:,idx_real);
        end

        pMSE_batch(ii) = mean((p_hat-c).^2);
    end

    pMSE_rep(reps) = pMSE_batch;

    % Kill the workers immediately after this small batch
    delete(gcp('nocreate'));

    fprintf('  Completed repetitions %d-%d\n',batch_start,batch_end);
end

mean_pMSE = mean(pMSE_rep);
sd_pMSE = std(pMSE_rep);

end