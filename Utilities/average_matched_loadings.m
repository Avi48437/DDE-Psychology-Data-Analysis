function out = average_matched_loadings(results, tau)

%% Default threshold

if nargin < 2
    tau = 0;
end

%% Identify completed iterations

valid_iteration = find(arrayfun( ...
    @(x) ~isempty(x.B1_CSP) && ~isempty(x.B2_CSP), results));

if isempty(valid_iteration)
    error('No completed DDE fits were found in results.');
end

reference_iteration = valid_iteration(1);

B1_reference = results(reference_iteration).B1_CSP;
B2_reference = results(reference_iteration).B2_CSP;

J = size(B1_reference, 1);
K1 = size(B1_reference, 2) - 1;
K2 = size(B2_reference, 2) - 1;

number_of_valid_iterations = length(valid_iteration);

if size(B2_reference, 1) ~= K1
    error('The number of B2 rows must equal the number of B1 factors.');
end

%% Storage

matched_B1 = NaN(J, K1 + 1, number_of_valid_iterations);
matched_B2 = NaN(K1, K2 + 1, number_of_valid_iterations);

B1_assignments = NaN(number_of_valid_iterations, K1);
B2_assignments = NaN(number_of_valid_iterations, K2);

B1_signs = ones(number_of_valid_iterations, K1);
B2_signs = ones(number_of_valid_iterations, K2);

%% Store reference fit

matched_B1(:, :, 1) = B1_reference;
matched_B2(:, :, 1) = B2_reference;

B1_assignments(1, :) = 1:K1;
B2_assignments(1, :) = 1:K2;

B1_reference_loadings = B1_reference(:, 2:end);
B2_reference_loadings = B2_reference(:, 2:end);

%% Match each remaining fit to iteration 1

for r = 2:number_of_valid_iterations

    iteration_number = valid_iteration(r);

    B1_current = results(iteration_number).B1_CSP;
    B2_current = results(iteration_number).B2_CSP;

    if ~isequal(size(B1_current), size(B1_reference))
        error('B1 size differs in iteration %d.', iteration_number);
    end

    if ~isequal(size(B2_current), size(B2_reference))
        error('B2 size differs in iteration %d.', iteration_number);
    end

    B1_current_intercept = B1_current(:, 1);
    B1_current_loadings = B1_current(:, 2:end);

    B2_current_intercept = B2_current(:, 1);
    B2_current_loadings = B2_current(:, 2:end);

    %% Match B1 columns

    B1_cost = zeros(K1, K1);

    for k_reference = 1:K1
        for k_current = 1:K1

            reference_column = ...
                B1_reference_loadings(:, k_reference);

            current_column = ...
                B1_current_loadings(:, k_current);

            denominator = ...
                norm(reference_column) * norm(current_column);

            if denominator == 0
                similarity = 0;
            else
                similarity = abs( ...
                    reference_column' * current_column / denominator);
            end

            B1_cost(k_reference, k_current) = 1 - similarity;
        end
    end

    B1_assignment = munkres(B1_cost);

    if any(B1_assignment == 0)
        error('Incomplete B1 assignment in iteration %d.', ...
            iteration_number);
    end

    B1_assignments(r, :) = B1_assignment;
    B1_reordered = B1_current_loadings(:, B1_assignment);

    %% Align B1 signs

    for k = 1:K1

        if B1_reference_loadings(:, k)' * B1_reordered(:, k) < 0
            B1_reordered(:, k) = -B1_reordered(:, k);
            B1_signs(r, k) = -1;
        end
    end

    matched_B1(:, :, r) = ...
        [B1_current_intercept, B1_reordered];

    %% Reorder B2 rows using the B1 assignment

    B2_intercept_reordered = ...
        B2_current_intercept(B1_assignment, :);

    B2_rows_reordered = ...
        B2_current_loadings(B1_assignment, :);

    %% Match B2 columns

    B2_cost = zeros(K2, K2);

    for k_reference = 1:K2
        for k_current = 1:K2

            reference_column = ...
                B2_reference_loadings(:, k_reference);

            current_column = ...
                B2_rows_reordered(:, k_current);

            denominator = ...
                norm(reference_column) * norm(current_column);

            if denominator == 0
                similarity = 0;
            else
                similarity = abs( ...
                    reference_column' * current_column / denominator);
            end

            B2_cost(k_reference, k_current) = 1 - similarity;
        end
    end

    B2_assignment = munkres(B2_cost);

    if any(B2_assignment == 0)
        error('Incomplete B2 assignment in iteration %d.', ...
            iteration_number);
    end

    B2_assignments(r, :) = B2_assignment;
    B2_reordered = B2_rows_reordered(:, B2_assignment);

    %% Align B2 signs

    for k = 1:K2

        if B2_reference_loadings(:, k)' * B2_reordered(:, k) < 0
            B2_reordered(:, k) = -B2_reordered(:, k);
            B2_signs(r, k) = -1;
        end
    end

    matched_B2(:, :, r) = ...
        [B2_intercept_reordered, B2_reordered];
end

%% Average matched matrices

average_B1_raw = mean(matched_B1, 3, 'omitnan');
average_B2_raw = mean(matched_B2, 3, 'omitnan');

%% Apply tau to non-intercept loading entries

average_B1_tau = average_B1_raw;
average_B2_tau = average_B2_raw;

average_B1_tau(:, 2:end) = ...
    thres(average_B1_tau(:, 2:end), tau);

average_B2_tau(:, 2:end) = ...
    thres(average_B2_tau(:, 2:end), tau);

%% Zero B2 rows corresponding to zero B1 columns

inactive_B1_columns = ...
    all(average_B1_tau(:, 2:end) == 0, 1);

average_B2_tau(inactive_B1_columns, :) = 0;

%% Identify active dimensions after thresholding

active_B1_columns = find(~inactive_B1_columns);

inactive_B2_columns = ...
    all(average_B2_tau(:, 2:end) == 0, 1);

active_B2_columns = find(~inactive_B2_columns);

%% Return results

out.reference_iteration = reference_iteration;
out.valid_iterations = valid_iteration;
out.tau = tau;

out.average_B1_raw = average_B1_raw;
out.average_B2_raw = average_B2_raw;

out.average_B1 = average_B1_tau;
out.average_B2 = average_B2_tau;

out.active_B1_columns = active_B1_columns;
out.inactive_B1_columns = find(inactive_B1_columns);

out.active_B2_columns = active_B2_columns;
out.inactive_B2_columns = find(inactive_B2_columns);

out.matched_B1 = matched_B1;
out.matched_B2 = matched_B2;

out.B1_assignments = B1_assignments;
out.B2_assignments = B2_assignments;

out.B1_signs = B1_signs;
out.B2_signs = B2_signs;

fprintf('Reference iteration: %d\n', reference_iteration);
fprintf('Number of averaged iterations: %d\n', ...
    number_of_valid_iterations);
fprintf('Applied threshold tau: %.4f\n', tau);
fprintf('Active average B1 columns: %d of %d\n', ...
    length(active_B1_columns), K1);
fprintf('Active average B2 columns: %d of %d\n', ...
    length(active_B2_columns), K2);

end