function [X_complete,Z_complete] = draw_completed_IPIP(X_obs,missing_mask,Z,A1_sample_long,B1_CSP,gamma_CSP,temp,z_cut)

[N,J] = size(X_obs);
C = size(A1_sample_long,3);

A1_draw = A1_sample_long(:,:,randi(C));
mu = A1_draw*B1_CSP.';

X_complete = X_obs;
Z_complete = Z;

for j = 1:J
    idx = missing_mask(:,j);
    n_mis = sum(idx);

    if n_mis == 0
        continue
    end

    sd_j = sqrt(gamma_CSP(j)/temp);
    Z_complete(idx,j) = mu(idx,j) + sd_j*randn(n_mis,1);
    X_complete(idx,j) = 1 + sum(Z_complete(idx,j) > z_cut(j,:),2);
end

if any(isnan(X_complete),'all')
    error('The candidate completed dataset still contains missing values.');
end

if any(~ismember(X_complete,1:5),'all')
    error('The candidate completed dataset contains values outside 1–5.');
end

end