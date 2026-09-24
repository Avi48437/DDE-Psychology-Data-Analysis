function Fit_out = threshold_model(Fit,threshold)

% THRESHOLD_MODEL
%
% Remove weak Layer-1 factors from a post-processed DDE fit using
% the relative L1 norm of each B1 loading column.
%
% For Layer-1 factor k:
%
%   L1_k = sum_j |B1(j,k)|
%
% and
%
%   RelativeL1_k = L1_k / max_l L1_l.
%
% Factor k is retained if:
%
%   RelativeL1_k >= threshold.
%
% The corresponding structures are pruned consistently:
%
%   B1         -> columns
%   B2         -> rows
%   A1         -> columns, if present
%   factor_map -> rows
%
% After Layer-1 pruning, any Layer-2 factor having an entirely
% zero loading column is also removed from:
%
%   B2
%   A2         -> columns, if present
%   A2_map     -> rows, if present
%
% Canonical factor numbering is rebuilt after pruning.
%
% Nothing is saved to disk.
%
% INPUTS
% ------
% Fit       : post-processed DDE fit
% threshold : relative L1 threshold in [0,1]
%
% OUTPUT
% -------
% Fit_out   : thresholded post-processed DDE fit
%
% EXAMPLE
% -------
% IPIPFMM = threshold_model(IPIPFMM_CC1,0.20);

%% ============================================================
% Validate
%% ============================================================

if nargin < 2
    error('threshold_model requires a fit and threshold.');
end

if ~isscalar(threshold) || threshold<0 || threshold>1
    error('threshold must be a scalar between 0 and 1.');
end

required_fields = {'B1','B2','factor_map'};

for j = 1:numel(required_fields)
    if ~isfield(Fit,required_fields{j})
        error('Fit is missing field %s.',required_fields{j});
    end
end

K1 = size(Fit.B1,2)-1;

if size(Fit.B2,1) ~= K1
    error('B2 rows do not match the number of Layer-1 factors.');
end

if height(Fit.factor_map) ~= K1
    error('factor_map rows do not match the number of Layer-1 factors.');
end

%% ============================================================
% Layer-1 L1 strengths
%% ============================================================

B1_loadings = Fit.B1(:,2:end);

l1_strength = sum(abs(B1_loadings),1);

max_l1 = max(l1_strength);

if max_l1==0
    error('All Layer-1 loading columns have zero L1 norm.');
end

relative_l1 = l1_strength/max_l1;

keep_A1 = relative_l1 >= threshold;

if ~any(keep_A1)
    error('Threshold %.3f removes every Layer-1 factor.',threshold);
end

old_K1 = K1;
old_canonical_A1 = find(keep_A1);

%% ============================================================
% Start output object
%% ============================================================

Fit_out = Fit;

%% ============================================================
% Prune B1
%% ============================================================

Fit_out.B1 = [ ...
    Fit.B1(:,1), ...
    Fit.B1(:,find(keep_A1)+1)];

%% ============================================================
% Prune B2 rows
%% ============================================================

Fit_out.B2 = Fit.B2(keep_A1,:);

%% ============================================================
% Prune A1 if available
%% ============================================================

if isfield(Fit,'A1') && ~isempty(Fit.A1)

    if size(Fit.A1,2) ~= old_K1
        error('A1 columns do not match the number of Layer-1 factors.');
    end

    Fit_out.A1 = Fit.A1(:,keep_A1);

end

%% ============================================================
% Prune factor map
%% ============================================================

factor_map = Fit.factor_map(keep_A1,:);

%% ============================================================
% Preserve old canonical identities
%% ============================================================

if ~ismember('PreThresholdCanonicalFactor', ...
        factor_map.Properties.VariableNames)

    factor_map.PreThresholdCanonicalFactor = ...
        old_canonical_A1(:);

else

    factor_map.PreThresholdCanonicalFactor = ...
        factor_map.CanonicalFactor;

end

%% ============================================================
% Rebuild canonical Layer-1 numbering
%% ============================================================

new_K1 = sum(keep_A1);

new_canonical = (1:new_K1)';

if ismember('CanonicalFactor',factor_map.Properties.VariableNames)
    factor_map.CanonicalFactor = new_canonical;
end

if ismember('CanonicalName',factor_map.Properties.VariableNames)
    factor_map.CanonicalName = "A1_" + string(new_canonical);
end

if ismember('Factor',factor_map.Properties.VariableNames)
    factor_map.Factor = "A1_" + string(new_canonical);
end

Fit_out.factor_map = factor_map;
Fit_out.K1 = new_K1;

%% ============================================================
% Recompute domain counts
%% ============================================================

if isfield(Fit_out,'domain_names')

    if ismember('Domain',factor_map.Properties.VariableNames)

        factor_domains = string(factor_map.Domain);

    elseif ismember('FactorDomain',factor_map.Properties.VariableNames)

        factor_domains = string(factor_map.FactorDomain);

    else

        factor_domains = strings(new_K1,1);

    end

    domain_counts = zeros(1,numel(Fit_out.domain_names));

    for d = 1:numel(Fit_out.domain_names)
        domain_counts(d) = sum( ...
            factor_domains==string(Fit_out.domain_names(d)));
    end

    Fit_out.domain_counts = domain_counts;

end

%% ============================================================
% Layer-2 pruning
%
% Remove any Layer-2 factor whose loading column becomes
% entirely zero after Layer-1 pruning.
%% ============================================================

old_K2 = size(Fit_out.B2,2)-1;

if old_K2 > 0

    B2_loadings = Fit_out.B2(:,2:end);

    keep_A2 = any(B2_loadings~=0,1);

    Fit_out.B2 = [ ...
        Fit_out.B2(:,1), ...
        B2_loadings(:,keep_A2)];

    new_K2 = sum(keep_A2);

    %% --------------------------------------------------------
    % Prune A2 if available
    %% --------------------------------------------------------

    if isfield(Fit,'A2') && ~isempty(Fit.A2)

        if size(Fit.A2,2) ~= old_K2
            error('A2 columns do not match the number of Layer-2 factors.');
        end

        Fit_out.A2 = Fit.A2(:,keep_A2);

    end

    %% --------------------------------------------------------
    % Prune A2 map if available
    %% --------------------------------------------------------

    if isfield(Fit,'A2_map') && ~isempty(Fit.A2_map)

        if height(Fit.A2_map) ~= old_K2
            error('A2_map rows do not match the number of Layer-2 factors.');
        end

        A2_map = Fit.A2_map(keep_A2,:);

        old_canonical_A2 = find(keep_A2);

        if ~ismember('PreThresholdCanonicalA2', ...
                A2_map.Properties.VariableNames)

            A2_map.PreThresholdCanonicalA2 = ...
                old_canonical_A2(:);

        end

        new_canonical_A2 = (1:new_K2)';

        if ismember('CanonicalA2',A2_map.Properties.VariableNames)
            A2_map.CanonicalA2 = new_canonical_A2;
        end

        if ismember('A2Factor',A2_map.Properties.VariableNames)
            A2_map.A2Factor = "A2_" + string(new_canonical_A2);
        end

        Fit_out.A2_map = A2_map;

    end

    Fit_out.K2 = new_K2;

else

    keep_A2 = false(1,0);
    new_K2 = 0;
    Fit_out.K2 = 0;

end

%% ============================================================
% Thresholding information
%% ============================================================

Fit_out.thresholding = struct;

Fit_out.thresholding.method = "Relative Layer-1 column L1 norm";
Fit_out.thresholding.threshold = threshold;

Fit_out.thresholding.original_K1 = old_K1;
Fit_out.thresholding.final_K1 = new_K1;

Fit_out.thresholding.L1Strength = l1_strength;
Fit_out.thresholding.RelativeL1 = relative_l1;

Fit_out.thresholding.KeepA1 = keep_A1;
Fit_out.thresholding.RemovedA1 = ~keep_A1;
Fit_out.thresholding.PreThresholdCanonicalA1 = old_canonical_A1;

Fit_out.thresholding.original_K2 = old_K2;
Fit_out.thresholding.final_K2 = new_K2;
Fit_out.thresholding.KeepA2 = keep_A2;

%% ============================================================
% Print summary
%% ============================================================

fprintf('\n============================================================\n');
fprintf('DDE MODEL THRESHOLDING\n');
fprintf('============================================================\n');
fprintf('Relative L1 threshold : %.3f\n',threshold);
fprintf('Layer 1: %d -> %d factors\n',old_K1,new_K1);
fprintf('Layer 2: %d -> %d factors\n',old_K2,new_K2);
fprintf('Retained old A1 indices: [%s]\n',num2str(old_canonical_A1));
fprintf('Removed old A1 indices : [%s]\n',num2str(find(~keep_A1)));

if isfield(Fit_out,'domain_counts')
    fprintf('Domain counts: [%s]\n',num2str(Fit_out.domain_counts));
end

fprintf('============================================================\n');

end