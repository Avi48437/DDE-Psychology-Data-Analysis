function [B1_aligned,B2_aligned,info] = ...
    align_combination_B1_B2(B1,B2,K1,K2,tol)

n_runs = size(B1,3);
J = size(B1,1);

B1_aligned = zeros(J,K1+1,n_runs);
B2_aligned = zeros(K1,K2+1,n_runs);

B1_assignments = cell(n_runs,1);
B2_assignments = cell(n_runs,1);


%% Reference repetition

B1_ref_full = B1(:,:,1);
B2_ref_full = B2(:,:,1);

active_B1_ref = find(any(abs(B1_ref_full(:,2:end)) > tol,1));

if numel(active_B1_ref) ~= K1
    error('Reference run has %d active B1 factors; expected %d.', ...
        numel(active_B1_ref),K1);
end

B1_ref = B1_ref_full(:,active_B1_ref+1);
B2_ref_rows = B2_ref_full(active_B1_ref,:);

if K2 > 0
    active_B2_ref = find(any(abs(B2_ref_rows(:,2:end)) > tol,1));

    if numel(active_B2_ref) ~= K2
        error('Reference run has %d active B2 factors; expected %d.', ...
            numel(active_B2_ref),K2);
    end
else
    active_B2_ref = [];
end

B1_aligned(:,:,1) = [B1_ref_full(:,1),B1_ref];

if K2 > 0
    B2_aligned(:,:,1) = ...
        [B2_ref_rows(:,1),B2_ref_rows(:,active_B2_ref+1)];
else
    B2_aligned(:,:,1) = B2_ref_rows(:,1);
end

B1_assignments{1} = active_B1_ref;
B2_assignments{1} = active_B2_ref;


%% Align the remaining repetitions

for r = 2:n_runs

    B1_r = B1(:,:,r);
    B2_r = B2(:,:,r);

    active_B1_r = find(any(abs(B1_r(:,2:end)) > tol,1));

    if numel(active_B1_r) ~= K1
        error('Run %d has %d active B1 factors; expected %d.', ...
            r,numel(active_B1_r),K1);
    end

    current_B1 = B1_r(:,active_B1_r+1);

    similarity_B1 = cosine_similarity_matrix(B1_ref,current_B1);

    assignment_B1 = munkres(-abs(similarity_B1));

    if any(assignment_B1 == 0)
        error('B1 matching failed in run %d.',r);
    end

    ordered_B1_ids = active_B1_r(assignment_B1);


    %% Align the binary-factor orientation exactly

    for k = 1:K1

        current_id = ordered_B1_ids(k);
        signed_similarity = ...
            similarity_B1(k,assignment_B1(k));

        if signed_similarity < 0

            loading = B1_r(:,current_id+1);

            % A1_new = 1 - A1_old:
            % b0_new = b0_old + b_old
            % b_new  = -b_old
            B1_r(:,1) = B1_r(:,1) + loading;
            B1_r(:,current_id+1) = -loading;

            % logit P(1-A1=1|A2) = -logit P(A1=1|A2)
            B2_r(current_id,:) = -B2_r(current_id,:);
        end
    end

    B1_aligned(:,:,r) = ...
        [B1_r(:,1),B1_r(:,ordered_B1_ids+1)];

    B2_r_ordered = B2_r(ordered_B1_ids,:);

    B1_assignments{r} = ordered_B1_ids;


    %% Align Layer-2 columns

    if K2 == 0
        B2_aligned(:,:,r) = B2_r_ordered(:,1);
        B2_assignments{r} = [];
        continue
    end

    active_B2_r = find( ...
        any(abs(B2_r_ordered(:,2:end)) > tol,1));

    if numel(active_B2_r) ~= K2
        error('Run %d has %d active B2 factors; expected %d.', ...
            r,numel(active_B2_r),K2);
    end

    reference_B2 = ...
        B2_ref_rows(:,active_B2_ref+1);

    current_B2 = ...
        B2_r_ordered(:,active_B2_r+1);

    similarity_B2 = ...
        cosine_similarity_matrix(reference_B2,current_B2);

    assignment_B2 = munkres(-abs(similarity_B2));

    if any(assignment_B2 == 0)
        error('B2 matching failed in run %d.',r);
    end

    ordered_B2_ids = active_B2_r(assignment_B2);


    %% Align orientation of each Layer-2 binary factor

    for h = 1:K2

        current_id = ordered_B2_ids(h);
        signed_similarity = ...
            similarity_B2(h,assignment_B2(h));

        if signed_similarity < 0

            loading = B2_r_ordered(:,current_id+1);

            % A2_new = 1 - A2_old
            B2_r_ordered(:,1) = ...
                B2_r_ordered(:,1) + loading;

            B2_r_ordered(:,current_id+1) = -loading;
        end
    end

    B2_aligned(:,:,r) = ...
        [B2_r_ordered(:,1), ...
         B2_r_ordered(:,ordered_B2_ids+1)];

    B2_assignments{r} = ordered_B2_ids;
end


%% Output alignment information

info.ReferenceRun = 1;
info.B1ReferenceIDs = active_B1_ref(:);
info.B2ReferenceIDs = active_B2_ref(:);
info.B1Assignments = B1_assignments;
info.B2Assignments = B2_assignments;

end


function similarity = cosine_similarity_matrix(A,B)

norm_A = sqrt(sum(A.^2,1));
norm_B = sqrt(sum(B.^2,1));

denominator = norm_A' * norm_B;
denominator(denominator == 0) = eps;

similarity = (A' * B) ./ denominator;

end