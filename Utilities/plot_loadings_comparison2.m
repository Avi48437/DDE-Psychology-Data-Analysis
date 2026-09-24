function plot_loadings_comparison2( ...
    B1_A,B2_A,groups_A,block_sizes_A,name_A, ...
    B1_B,B2_B,groups_B,block_sizes_B,name_B, ...
    trait_labels,order_by_B2)

%% ============================================================
% Defaults
% ============================================================

if nargin < 11 || isempty(trait_labels)
    trait_labels = {'EXT','EST','AGR','CSN','OPN'};
end

if nargin < 12 || isempty(order_by_B2)
    order_by_B2 = false;
end

zero_tol = 1e-12;

%% ============================================================
% Prepare both models
% ============================================================

[B1_plot_A,B2_plot_A,order_A,active2_A] = ...
    prepare_model(B1_A,B2_A,groups_A,block_sizes_A,zero_tol,order_by_B2);

[B1_plot_B,B2_plot_B,order_B,active2_B] = ...
    prepare_model(B1_B,B2_B,groups_B,block_sizes_B,zero_tol,order_by_B2);

%% ============================================================
% Common symmetric color scale
% ============================================================

max_abs_loading = max(abs([ ...
    B1_plot_A(:); ...
    B2_plot_A(:); ...
    B1_plot_B(:); ...
    B2_plot_B(:)]));

if max_abs_loading == 0
    max_abs_loading = 1;
end

color_limits = [-max_abs_loading max_abs_loading];

%% ============================================================
% Red-white-blue colormap
% ============================================================

number_of_colors = 256;
half_colors = floor(number_of_colors/2);

blue_to_white = [ ...
    linspace(0,1,half_colors)', ...
    linspace(0,1,half_colors)', ...
    ones(half_colors,1)];

white_to_red = [ ...
    ones(number_of_colors-half_colors,1), ...
    linspace(1,0,number_of_colors-half_colors)', ...
    linspace(1,0,number_of_colors-half_colors)'];

red_white_blue = [blue_to_white; white_to_red];

%% ============================================================
% Trait-block centers and boundaries
% ============================================================

edges_A = cumsum(block_sizes_A);
starts_A = [1 edges_A(1:end-1)+1];
centers_A = (starts_A + edges_A)/2;

edges_B = cumsum(block_sizes_B);
starts_B = [1 edges_B(1:end-1)+1];
centers_B = (starts_B + edges_B)/2;

%% ============================================================
% Plot
% ============================================================

figure('Position',[50 50 1600 1050],'Color','w');

tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

%% Model A: B1

nexttile

imagesc(B1_plot_A)
set(gca,'YDir','reverse')
axis tight
caxis(color_limits)
colormap(red_white_blue)
colorbar

title(sprintf('%s: B_1',name_A))
xlabel('Intercept and Layer 1 Factors')
ylabel('Personality Trait')

xticks(1:size(B1_plot_A,2))
xticklabels(["Intercept","A1_" + string(order_A)])
xtickangle(45)

yticks(centers_A)
yticklabels(trait_labels)

for j = 1:numel(edges_A)-1
    yline(edges_A(j)+0.5,'k-','LineWidth',1.2)
end

%% Model A: B2

nexttile

imagesc(B2_plot_A)
set(gca,'YDir','reverse')
axis tight
caxis(color_limits)
colormap(red_white_blue)
colorbar

title(sprintf('%s: B_2',name_A))
xlabel('Intercept and Layer 2 Factors')
ylabel('Layer 1 Factors')

xticks(1:size(B2_plot_A,2))
xticklabels(["Intercept","A2_" + string(active2_A)])

yticks(1:numel(order_A))
yticklabels("A1_" + string(order_A))

%% Model B: B1

nexttile

imagesc(B1_plot_B)
set(gca,'YDir','reverse')
axis tight
caxis(color_limits)
colormap(red_white_blue)
colorbar

title(sprintf('%s: B_1',name_B))
xlabel('Intercept and Layer 1 Factors')
ylabel('Personality Trait')

xticks(1:size(B1_plot_B,2))
xticklabels(["Intercept","A1_" + string(order_B)])
xtickangle(45)

yticks(centers_B)
yticklabels(trait_labels)

for j = 1:numel(edges_B)-1
    yline(edges_B(j)+0.5,'k-','LineWidth',1.2)
end

%% Model B: B2

nexttile

imagesc(B2_plot_B)
set(gca,'YDir','reverse')
axis tight
caxis(color_limits)
colormap(red_white_blue)
colorbar

title(sprintf('%s: B_2',name_B))
xlabel('Intercept and Layer 2 Factors')
ylabel('Layer 1 Factors')

xticks(1:size(B2_plot_B,2))
xticklabels(["Intercept","A2_" + string(active2_B)])

yticks(1:numel(order_B))
yticklabels("A1_" + string(order_B))

%% Overall title

if order_by_B2
    title(tl,'DDE Loading Matrix Comparison — Layer 1 ordered by B_2 sign')
else
    title(tl,'DDE Loading Matrix Comparison')
end

end


%% ============================================================
% Helper function
% ============================================================

function [B1_plot,B2_plot,order,active2] = ...
    prepare_model(B1,B2,groups,block_sizes,zero_tol,order_by_B2)

%% ============================================================
% Manual Layer-1 grouping/order
% ============================================================

manual_order = [groups{:}];

K1_max = size(B1,2)-1;

if any(manual_order < 1) || any(manual_order > K1_max)
    error('A supplied Layer-1 factor exceeds K1_max.');
end

if numel(unique(manual_order)) ~= numel(manual_order)
    error('A Layer-1 factor appears more than once in groups.');
end

%% ============================================================
% Check omitted Layer-1 factors
% ============================================================

all_factors = 1:K1_max;
omitted = setdiff(all_factors,manual_order);

if ~isempty(omitted)

    omitted_nonzero = omitted(any(abs(B1(:,omitted+1)) > zero_tol,1));

    if ~isempty(omitted_nonzero)
        error('Unlisted Layer-1 factors are nonzero: %s', ...
            mat2str(omitted_nonzero));
    end

end

%% ============================================================
% Check row blocks
% ============================================================

if numel(groups) ~= numel(block_sizes)
    error('Number of factor groups must equal number of row blocks.');
end

if sum(block_sizes) ~= size(B1,1)
    error('block_sizes must sum to the number of B1 rows.');
end

%% ============================================================
% Reorder observed-variable rows within each trait block
%
% First factor in each supplied group is the anchor.
%
% Rows:
%   positive anchor loading
%   zero anchor loading
%   negative anchor loading
% ============================================================

block_end = cumsum(block_sizes);
block_start = [1 block_end(1:end-1)+1];

row_order = [];

for d = 1:numel(groups)

    rows = block_start(d):block_end(d);
    anchor = groups{d}(1);

    loading = B1(rows,anchor+1);

    pos_idx = find(loading > zero_tol);
    zero_idx = find(abs(loading) <= zero_tol);
    neg_idx = find(loading < -zero_tol);

    % Positive: largest positive first
    [~,ii] = sort(loading(pos_idx),'descend');
    pos_idx = pos_idx(ii);

    % Negative: most negative first
    [~,ii] = sort(loading(neg_idx),'ascend');
    neg_idx = neg_idx(ii);

    block_order = [pos_idx; zero_idx; neg_idx];

    row_order = [row_order rows(block_order)];

end

%% ============================================================
% Find active Layer-2 factors
% ============================================================

B2_manual = B2(manual_order,:);

active2 = find(any(abs(B2_manual(:,2:end)) > zero_tol,1));

%% ============================================================
% Optional Layer-1 ordering according to B2 sign
% ============================================================

order = manual_order;

if order_by_B2

    if isempty(active2)

        warning('No active Layer-2 factor found. Keeping manual Layer-1 ordering.');

    elseif numel(active2) > 1

        error(['B2-sign ordering requires exactly one active Layer-2 ', ...
            'factor. This model has %d active Layer-2 factors.'], ...
            numel(active2));

    else

        % B2 column corresponding to the single active Layer-2 factor
        b2_factor = active2(1);

        % Values in the CURRENT manual Layer-1 order
        b2_loading = B2(manual_order,b2_factor+1);

        positive_idx = find(b2_loading > zero_tol);
        negative_idx = find(b2_loading < -zero_tol);
        zero_idx = find(abs(b2_loading) <= zero_tol);

        % Preserve manual ordering within each sign group
        sign_permutation = [positive_idx; negative_idx; zero_idx];

        order = manual_order(sign_permutation);

    end

end

%% ============================================================
% Construct B1
% Intercept retained
% ============================================================

B1_plot = [ ...
    B1(row_order,1), ...
    B1(row_order,order+1)];

%% ============================================================
% Construct B2
% Same Layer-1 permutation as B1 columns
% ============================================================

B2_rows = B2(order,:);

% Recalculate active Layer-2 factors after row permutation
active2 = find(any(abs(B2_rows(:,2:end)) > zero_tol,1));

B2_plot = [ ...
    B2_rows(:,1), ...
    B2_rows(:,active2+1)];

end