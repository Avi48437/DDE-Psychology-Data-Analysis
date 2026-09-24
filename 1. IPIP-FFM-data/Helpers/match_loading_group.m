function matched_group = match_loading_group(results_subset, repetition_ids)

% MATCH_LOADING_GROUP
%
% Matches all B1 and B2 loading matrices in one fixed (K1,K2)
% combination to the first repetition in that group.
%
% OUTPUT
% ------
% matched_group.B1
%     J x (K1+1) x number_of_runs
%
% matched_group.B2
%     K1 x (K2+1) x number_of_runs
%
% The first column remains the intercept.
% Active factors are placed in the same order as the reference run.

n_runs = numel(results_subset);

if nargin < 2
    repetition_ids = (1:n_runs)';
end

%% Reference run

reference_B1 = results_subset(1).B1_CSP;
reference_B2 = results_subset(1).B2_CSP;

reference_B1_active = ...
    find(any(reference_B1(:, 2:end) ~= 0, 1)) + 1;

K1 = numel(reference_B1_active);

reference_B2_rows = reference_B1_active - 1;

reference_B2_active = ...
    find(any(reference_B2(reference_B2_rows, 2:end) ~= 0, 1)) + 1;

K2 = numel(reference_B2_active);

%% Compact reference matrices

reference_B1_compact = [ ...
    reference_B1(:, 1), ...
    reference_B1(:, reference_B1_active)];

reference_B2_compact = [ ...
    reference_B2(reference_B2_rows, 1), ...
    reference_B2(reference_B2_rows, reference_B2_active)];

%% Storage

J = size(reference_B1, 1);

matched_B1 = zeros(J, K1 + 1, n_runs);
matched_B2 = zeros(K1, K2 + 1, n_runs);

B1_permutation = zeros(n_runs, K1);
B2_permutation = zeros(n_runs, K2);

matched_B1(:, :, 1) = reference_B1_compact;
matched_B2(:, :, 1) = reference_B2_compact;

B1_permutation(1, :) = reference_B1_active - 1;

if K2 > 0
    B2_permutation(1, :) = reference_B2_active - 1;
end

%% Match every remaining run to the reference

for s = 2:n_runs

    B1 = results_subset(s).B1_CSP;
    B2 = results_subset(s).B2_CSP;

    active_B1 = find(any(B1(:, 2:end) ~= 0, 1)) + 1;

    if numel(active_B1) ~= K1
        error( ...
            'Run %d has %d active B1 factors; expected %d.', ...
            repetition_ids(s), numel(active_B1), K1);
    end

    %% Match B1 columns

    reference_loadings = reference_B1(:, reference_B1_active);
    current_loadings = B1(:, active_B1);

    similarity_B1 = cosine_similarity_matrix( ...
        reference_loadings, current_loadings);

    pairs_B1 = matchpairs(-similarity_B1, 1e10);
    pairs_B1 = sortrows(pairs_B1, 1);

    B1_order = pairs_B1(:, 2)';
    matched_active_B1 = active_B1(B1_order);

    matched_B1(:, :, s) = [ ...
        B1(:, 1), ...
        B1(:, matched_active_B1)];

    B1_permutation(s, :) = matched_active_B1 - 1;

    %% Reorder B2 rows according to the matched B1 factors

    matched_B2_rows = matched_active_B1 - 1;
    B2_row_matched = B2(matched_B2_rows, :);

    active_B2 = ...
        find(any(B2_row_matched(:, 2:end) ~= 0, 1)) + 1;

    if numel(active_B2) ~= K2
        error( ...
            'Run %d has %d active B2 factors; expected %d.', ...
            repetition_ids(s), numel(active_B2), K2);
    end

    %% Match B2 columns

    if K2 == 0

        matched_B2(:, :, s) = B2_row_matched(:, 1);

    else

        reference_B2_loadings = ...
            reference_B2_compact(:, 2:end);

        current_B2_loadings = ...
            B2_row_matched(:, active_B2);

        similarity_B2 = cosine_similarity_matrix( ...
            reference_B2_loadings, current_B2_loadings);

        pairs_B2 = matchpairs(-similarity_B2, 1e10);
        pairs_B2 = sortrows(pairs_B2, 1);

        B2_order = pairs_B2(:, 2)';
        matched_active_B2 = active_B2(B2_order);

        matched_B2(:, :, s) = [ ...
            B2_row_matched(:, 1), ...
            B2_row_matched(:, matched_active_B2)];

        B2_permutation(s, :) = matched_active_B2 - 1;

    end

end

%% Output

matched_group.K1 = K1;
matched_group.K2 = K2;

matched_group.RepetitionIDs = repetition_ids(:);

matched_group.ReferenceRepetition = repetition_ids(1);

matched_group.B1 = matched_B1;
matched_group.B2 = matched_B2;

matched_group.B1Permutation = B1_permutation;
matched_group.B2Permutation = B2_permutation;

end


function similarity = cosine_similarity_matrix(reference, current)

n_reference = size(reference, 2);
n_current = size(current, 2);

similarity = zeros(n_reference, n_current);

for i = 1:n_reference

    x = reference(:, i);

    for j = 1:n_current

        y = current(:, j);

        similarity(i, j) = ...
            dot(x, y) / (norm(x) * norm(y) + eps);

    end

end

end