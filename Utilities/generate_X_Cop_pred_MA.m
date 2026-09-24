function [X_sim,Z_sim,margin] = generate_X_Cop_pred_MA(N,prop_true,B1,B2,gamma,X_data_obs,Z_train)
% GENERATE_X_COP_PRED_MA
% Synthetic generation using the Feldman-Kowal margin adjustment.
%
% For each item j:
%
%   Z_j^n(x) = max[{Z_ij^obs:X_ij^obs<=x}
%                  union
%                  {Z_ij^obs:X_ij^obs=min(X_j^obs)}]
%
%   F_tilde_j(x) = G_j{Z_j^n(x)}
%
% where G_j is the exact marginal CDF induced by the fitted DDE.
%
% IPIP support is assumed to be {1,2,3,4,5}.

    [n_train,p] = size(X_data_obs);
    prop_true = prop_true(:).';
    gamma = gamma(:).';
    K1 = size(B1,2)-1;
    K2 = numel(prop_true);
    x_support = (1:5)';
    n_levels = numel(x_support);

    if ~isequal(size(Z_train),[n_train,p]), error('Z_train must have the same size as X_data_obs.'); end
    if size(B1,1) ~= p, error('B1 must have one row for each item.'); end
    if size(B2,1) ~= K1 || size(B2,2) ~= K2+1, error('B2 has incompatible dimensions.'); end
    if numel(gamma) ~= p || any(~isfinite(gamma) | gamma <= 0), error('gamma must contain p positive finite values.'); end
    if any(~isfinite(prop_true) | prop_true < 0 | prop_true > 1), error('prop_true must contain probabilities in [0,1].'); end

    %% Generate synthetic latent data from the fitted DDE
    A2_sim = double(rand(N,K2) < prop_true);
    prob_A1_sim = logistic([ones(N,1),A2_sim]*B2');
    A1_sim = double(rand(N,K1) < prob_A1_sim);
    mu_sim = [ones(N,1),A1_sim]*B1';
    Z_sim = mu_sim+randn(N,p).*sqrt(gamma);

    %% Construct margin adjustment and transform Z_sim
    X_sim = zeros(N,p);

    margin.x_support = cell(1,p);
    margin.x_min_observed = zeros(1,p);
    margin.z_cut = cell(1,p);
    margin.F_equation9 = cell(1,p);
    margin.F_adjusted = cell(1,p);
    margin.category_prob = cell(1,p);
    margin.mixture_mean = cell(1,p);
    margin.mixture_weight = cell(1,p);

    for j = 1:p
        observed = ismember(X_data_obs(:,j),x_support) & isfinite(Z_train(:,j));
        x_obs = X_data_obs(observed,j);
        z_obs = Z_train(observed,j);

        if isempty(x_obs), error('Item %d contains no observed response in 1:5.',j); end

        %% Exact equation (10), including the union term
        x_min = min(x_obs);
        z_at_min = z_obs(x_obs == x_min);
        z_cut = zeros(n_levels,1);

        for r = 1:n_levels
            first_set = z_obs(x_obs <= x_support(r));
            second_set = z_at_min;
            z_cut(r) = max([first_set;second_set]);
        end

        if any(diff(z_cut) < 0), error('Equation (10) produced unordered cutpoints for item %d.',j); end

        %% Exact model-induced marginal CDF G_j
        [mixture_mean,mixture_weight] = dde_item_mixture(j,prop_true,B1,B2);
        F_equation9 = dde_marginal_cdf(z_cut,mixture_mean,mixture_weight,sqrt(gamma(j)));

        if any(diff(F_equation9) < -1e-10), error('Equation (9) produced a decreasing CDF for item %d.',j); end

        %% IPIP has bounded support {1,2,3,4,5}
        % The first four cumulative probabilities define the five intervals:
        %
        % category 1:        0 < U <= F(1)
        % category 2:     F(1) < U <= F(2)
        % category 3:     F(2) < U <= F(3)
        % category 4:     F(3) < U <= F(4)
        % category 5:     F(4) < U <= 1
        %
        % Therefore F(5)=1 is imposed by the known upper support boundary.

        F_adjusted = [F_equation9(1:n_levels-1);1];
        category_prob = [F_adjusted(1);diff(F_adjusted)];

        if any(category_prob < -1e-10), error('Negative category probability produced for item %d.',j); end

        category_prob = max(category_prob,0);
        category_prob = category_prob/sum(category_prob);

        %% Exact probability-integral transform U = G_j(Z_sim)
        U_sim = dde_marginal_cdf(Z_sim(:,j),mixture_mean,mixture_weight,sqrt(gamma(j)));

        %% Generalized inverse of the adjusted discrete CDF
        idx_sim = 1+sum(U_sim > F_adjusted(1:end-1).',2);

        if any(idx_sim < 1 | idx_sim > n_levels), error('Inverse CDF failed for item %d.',j); end

        X_sim(:,j) = x_support(idx_sim);

        margin.x_support{j} = x_support;
        margin.x_min_observed(j) = x_min;
        margin.z_cut{j} = z_cut;
        margin.F_equation9{j} = F_equation9;
        margin.F_adjusted{j} = F_adjusted;
        margin.category_prob{j} = category_prob;
        margin.mixture_mean{j} = mixture_mean;
        margin.mixture_weight{j} = mixture_weight;
    end
end

function [mixture_mean,mixture_weight] = dde_item_mixture(j,prop_true,B1,B2)
% Construct the exact finite Gaussian-mixture marginal for Z_j.

    active_A1 = find(B1(j,2:end) ~= 0);

    if isempty(active_A1)
        mixture_mean = B1(j,1);
        mixture_weight = 1;
        return
    end

    active_A2 = find(any(B2(active_A1,2:end) ~= 0,1));
    A1_states = binary_states(numel(active_A1));
    A2_states = binary_states(numel(active_A2));

    n_A1 = size(A1_states,1);
    n_A2 = size(A2_states,1);

    mixture_mean = zeros(n_A1*n_A2,1);
    mixture_weight = zeros(n_A1*n_A2,1);
    position = 1;

    for h = 1:n_A2
        a2 = A2_states(h,:);

        if isempty(active_A2)
            weight_A2 = 1;
            eta_A1 = B2(active_A1,1);
        else
            p2 = prop_true(active_A2);
            weight_A2 = prod(a2.*p2+(1-a2).*(1-p2));
            eta_A1 = B2(active_A1,1)+B2(active_A1,active_A2+1)*a2';
        end

        prob_A1 = logistic(eta_A1).';
        weight_A1 = prod(A1_states.*prob_A1+(1-A1_states).*(1-prob_A1),2);
        mean_A1 = B1(j,1)+A1_states*B1(j,active_A1+1).';

        idx = position:(position+n_A1-1);
        mixture_mean(idx) = mean_A1;
        mixture_weight(idx) = weight_A2*weight_A1;
        position = position+n_A1;
    end

    keep = mixture_weight > 0;
    mixture_mean = mixture_mean(keep);
    mixture_weight = mixture_weight(keep);
    mixture_weight = mixture_weight/sum(mixture_weight);
end

function F = dde_marginal_cdf(z,mixture_mean,mixture_weight,sigma)
% Evaluate G_j(z) exactly as a finite Gaussian-mixture CDF.

    original_size = size(z);
    z = z(:);
    mixture_mean = mixture_mean(:);
    mixture_weight = mixture_weight(:);

    F = zeros(numel(z),1);
    block_size = 256;

    for first = 1:block_size:numel(mixture_mean)
        last = min(first+block_size-1,numel(mixture_mean));
        idx = first:last;
        standardized = (z-mixture_mean(idx).')/sigma;
        F = F+normcdf(standardized)*mixture_weight(idx);
    end

    F = min(max(F,0),1);
    F = reshape(F,original_size);
end

function states = binary_states(K)
% Enumerate every binary vector of length K.

    if K == 0
        states = zeros(1,0);
        return
    end

    n_states = 2^K;
    codes = uint64((0:n_states-1)');
    states = zeros(n_states,K);

    for k = 1:K
        states(:,k) = bitget(codes,k);
    end
end

function y = logistic(x)
% Numerically stable logistic function.

    y = zeros(size(x));
    positive = x >= 0;
    y(positive) = 1./(1+exp(-x(positive)));
    ex = exp(x(~positive));
    y(~positive) = ex./(1+ex);
end