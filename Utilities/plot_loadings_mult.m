function plot_loadings_mult(dataset_names,add_cluster_panels)

% ============================================================
% Plot one, two, or three post-processed DDE fits
%
% INPUTS
% -------
% dataset_names       : string / char / cell array of names of the
%                       post-processed structs already loaded in workspace
%
% add_cluster_panels  : true  -> B1, B2, clustered B1, clustered B2
%                       false -> B1 and B2 only
%
% Examples:
%
% plot_loadings_mult('IPIP98',false)
% plot_loadings_mult({'IPIP98','LOOPR'},true)
% plot_loadings_mult({'IPIPFMM','IPIP98','LOOPR'},true)
% ============================================================

if nargin < 2
    add_cluster_panels = false;
end

if ischar(dataset_names) || isstring(dataset_names)
    dataset_names = cellstr(dataset_names);
end

n_sets = numel(dataset_names);

if n_sets < 1 || n_sets > 3
    error('dataset_names must contain 1, 2, or 3 dataset names.');
end

D = cell(n_sets,1);

for s = 1:n_sets
    D{s} = evalin('base',dataset_names{s});
end

%% ============================================================
% Common color range
%% ============================================================

all_vals = [];

for s = 1:n_sets
    all_vals = [all_vals;D{s}.B1(:);D{s}.B2(:)];
end

cmax = max(abs(all_vals));

if cmax==0
    cmax = 1;
end

%% ============================================================
% Layout width
%% ============================================================

if add_cluster_panels
    n_cols = 6;
else
    n_cols = 3;
end

%% ============================================================
% Figure
%
% IMPORTANT:
% Do not specify Position, OuterPosition, WindowState, or
% WindowStyle here. This allows MATLAB to place this figure
% inside the common Figure Container as a tab.
%% ============================================================

fig = figure( ...
    'Color','w', ...
    'Name','Loading Matrices', ...
    'NumberTitle','off');

T = tiledlayout( ...
    fig, ...
    n_sets,n_cols, ...
    'TileSpacing','compact', ...
    'Padding','compact');

%% ============================================================
% Plot datasets
%% ============================================================

for s = 1:n_sets

    ds = D{s};
    row0 = (s-1)*n_cols;

    %% --------------------------------------------------------
    % Original B1
    %% --------------------------------------------------------

    ax1 = nexttile(T,row0+1,[1 2]);
    plot_B1_panel(ax1,ds,false,cmax);

    %% --------------------------------------------------------
    % Original B2
    %% --------------------------------------------------------

    ax2 = nexttile(T,row0+3,[1 1]);
    plot_B2_panel(ax2,ds,false,cmax);

    %% --------------------------------------------------------
    % Clustered panels
    %% --------------------------------------------------------

    if add_cluster_panels

        ax3 = nexttile(T,row0+4,[1 2]);
        plot_B1_panel(ax3,ds,true,cmax);

        ax4 = nexttile(T,row0+6,[1 1]);
        plot_B2_panel(ax4,ds,true,cmax);

    end

end

%% ============================================================
% Force all ordinary text black
%% ============================================================

set(findall(fig,'Type','text'),'Color','k');

end


% ============================================================
% B1 panel
% ============================================================

function plot_B1_panel(ax,ds,use_clustered,cmax)

B1 = ds.B1;
label = get_dataset_label(ds);

if use_clustered

    [ord,pos_n,zero_n,neg_n] = get_cluster_order(ds);

    B1 = [B1(:,1),B1(:,ord+1)];

    xlabels = ['Int',string(ord)];

    group_bounds = cumsum([pos_n,zero_n,neg_n]);

    ttl = sprintf('%s: B_1 (clustered)',label);

else

    K1 = size(B1,2)-1;

    xlabels = ['Int',string(1:K1)];

    group_bounds = get_domain_group_bounds(ds);

    ttl = sprintf('%s: B_1',label);

end

imagesc(ax,B1);

colormap(ax,red_white_blue(256));

caxis(ax,[-cmax cmax]);

cb = colorbar(ax);
cb.Color = 'k';
cb.TickDirection = 'out';

set( ...
    ax, ...
    'Box','on', ...
    'LineWidth',1, ...
    'FontSize',12, ...
    'XColor','k', ...
    'YColor','k', ...
    'Layer','top', ...
    'TickDir','out', ...
    'Color','w');

title( ...
    ax, ...
    ttl, ...
    'FontSize',16, ...
    'FontWeight','bold', ...
    'Color','k');

xticks(ax,1:size(B1,2));
xticklabels(ax,xlabels);
xtickangle(ax,45);

[y_ticks,y_labels,y_bounds] = get_row_info(ds);

yticks(ax,y_ticks);
yticklabels(ax,y_labels);

ylabel( ...
    ax, ...
    'Personality Trait', ...
    'FontSize',15, ...
    'Color','k');

hold(ax,'on');

for b = y_bounds
    yline(ax,b+0.5,'k-','LineWidth',1);
end

xline(ax,1.5,'k-','LineWidth',1);

for b = group_bounds

    if b < size(B1,2)-1
        xline(ax,1.5+b,'k-','LineWidth',1);
    end

end

hold(ax,'off');

axis(ax,'tight');

end


% ============================================================
% B2 panel
% ============================================================

function plot_B2_panel(ax,ds,use_clustered,cmax)

B2 = ds.B2;
label = get_dataset_label(ds);

if use_clustered

    [ord,pos_n,zero_n,neg_n] = get_cluster_order(ds);

    B2 = B2(ord,:);

    ylabels = string(ord);

    row_bounds = cumsum([pos_n,zero_n,neg_n]);

    ttl = sprintf('%s: B_2 (clustered)',label);

else

    K1 = size(B2,1);

    ylabels = string(1:K1);

    row_bounds = get_domain_group_bounds(ds);

    ttl = sprintf('%s: B_2',label);

end

K2 = size(B2,2)-1;

xlabels = ['Int',string(1:K2)];

imagesc(ax,B2);

colormap(ax,red_white_blue(256));

caxis(ax,[-cmax cmax]);

cb = colorbar(ax);
cb.Color = 'k';
cb.TickDirection = 'out';

set( ...
    ax, ...
    'Box','on', ...
    'LineWidth',1, ...
    'FontSize',12, ...
    'XColor','k', ...
    'YColor','k', ...
    'Layer','top', ...
    'TickDir','out', ...
    'Color','w');

title( ...
    ax, ...
    ttl, ...
    'FontSize',16, ...
    'FontWeight','bold', ...
    'Color','k');

xticks(ax,1:size(B2,2));
xticklabels(ax,xlabels);
xtickangle(ax,45);

yticks(ax,1:size(B2,1));
yticklabels(ax,ylabels);

ylabel( ...
    ax, ...
    'Layer 1 Factors', ...
    'FontSize',15, ...
    'Color','k');

hold(ax,'on');

for b = row_bounds

    if b < size(B2,1)
        yline(ax,b+0.5,'k-','LineWidth',1);
    end

end

xline(ax,1.5,'k-','LineWidth',1);

hold(ax,'off');

axis(ax,'tight');

end


% ============================================================
% Row ticks and labels for B1
% ============================================================

function [y_ticks,y_labels,y_bounds] = get_row_info(ds)

if isfield(ds,'block_sizes')

    block_sizes = ds.block_sizes(:)';

elseif isfield(ds,'block_size') && isfield(ds,'domain_names')

    block_sizes = repmat( ...
        ds.block_size, ...
        1, ...
        numel(ds.domain_names));

else

    error('Could not find block_sizes information in dataset.');

end

if isfield(ds,'domain_names')

    y_labels = cellstr(ds.domain_names(:));

else

    y_labels = cellstr( ...
        "Block " + string(1:numel(block_sizes)));

end

y_ticks = zeros(1,numel(block_sizes));

start_idx = 1;

for j = 1:numel(block_sizes)

    y_ticks(j) = ...
        start_idx + ...
        (block_sizes(j)-1)/2;

    start_idx = ...
        start_idx + block_sizes(j);

end

y_bounds = cumsum(block_sizes);

y_bounds = y_bounds(1:end-1);

end


% ============================================================
% Domain group bounds for original factor ordering
% ============================================================

function bounds = get_domain_group_bounds(ds)

if isfield(ds,'domain_counts')

    counts = ds.domain_counts(:)';

    bounds = cumsum(counts);

    bounds = bounds(1:end-1);

else

    bounds = [];

end

end


% ============================================================
% Cluster ordering from B2
%
% positive -> zero -> negative
% ============================================================

function [ord,pos_n,zero_n,neg_n] = get_cluster_order(ds)

B2 = ds.B2;

if size(B2,2)<2

    ord = 1:size(B2,1);

    pos_n = size(B2,1);
    zero_n = 0;
    neg_n = 0;

    return

end

score = B2(:,2);

pos_idx = find(score>0);
zero_idx = find(score==0);
neg_idx = find(score<0);

[~,pord] = sort(score(pos_idx),'descend');
[~,zord] = sort(score(zero_idx),'descend');
[~,nord] = sort(score(neg_idx),'descend');

pos_idx = pos_idx(pord);
zero_idx = zero_idx(zord);
neg_idx = neg_idx(nord);

ord = [pos_idx;zero_idx;neg_idx]';

pos_n = numel(pos_idx);
zero_n = numel(zero_idx);
neg_n = numel(neg_idx);

end


% ============================================================
% Dataset label
% ============================================================

function label = get_dataset_label(ds)

if isfield(ds,'dataset_name')

    label = char(ds.dataset_name);

else

    label = 'Dataset';

end

end


% ============================================================
% Red-white-blue colormap
% ============================================================

function cmap = red_white_blue(m)

if nargin<1
    m = 256;
end

m1 = floor(m/2);
m2 = m-m1;

blue_to_white = [ ...
    linspace(0,1,m1)', ...
    linspace(0,1,m1)', ...
    ones(m1,1)];

white_to_red = [ ...
    ones(m2,1), ...
    linspace(1,0,m2)', ...
    linspace(1,0,m2)'];

cmap = [blue_to_white;white_to_red];

end