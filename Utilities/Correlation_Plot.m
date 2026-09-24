%% ============================================================
% Survey-question correlation matrices
%
% X = observed questionnaire responses used for DDE fitting
%
% Item ordering:
%   exactly the same ordering used in the final B1 loading plots
%
% IPIP-100 : one fixed dataset
% IPIP-FMM : average correlation across 100 fitted samples
% LOOPR    : one fixed dataset
%% ============================================================

clear
clc
close all

root_dir = '/Users/avinandanroy/Desktop/Research/DDE/DDE_Fitting_Pshycology';

addpath(fullfile(root_dir,'Scripts'));
addpath(fullfile(root_dir,'Algorithms'));
addpath(fullfile(root_dir,'Utilities'));
addpath(fullfile(root_dir,'1. IPIP-FFM-data'));

%% ============================================================
% 1. IPIP-100
%% ============================================================

load(fullfile(root_dir,'2. IPIP 100','Analysis', ...
    'IPIP100_PostProcessed.mat'),'IPIP100');

B5_wo = readtable( ...
    fullfile(root_dir,'2. IPIP 100','B5_wo.csv'), ...
    'VariableNamingRule','preserve');

big5_names = "BIGFIVE_" + string(1:100);

assert(all(ismember(big5_names,string(B5_wo.Properties.VariableNames))), ...
    'Some BIGFIVE_1,...,BIGFIVE_100 variables are missing.');

X100_original = B5_wo{:,big5_names};

% Exact final B1 item ordering
row_order_100 = IPIP100.item_map.OriginalRow(:)';

assert(isequal(row_order_100,IPIP100.row_order(:)'), ...
    'IPIP100 item_map and row_order disagree.');

X100 = X100_original(:,row_order_100);

R100 = corr(X100,'Rows','pairwise');

fprintf('\nIPIP-100 X size: %s\n',mat2str(size(X100)));

%% ============================================================
% 2. LOOPR
%
% LOOPR.X has already received the exact same item permutation
% as LOOPR.B1 during postprocessing.
%% ============================================================

load(fullfile(root_dir,'LOPR','Analysis', ...
    'LOOPR_PostProcessed.mat'),'LOOPR');

XLOOPR = double(LOOPR.X);

assert(size(XLOOPR,2)==60);
assert(size(XLOOPR,2)==size(LOOPR.B1,1));

RLOOPR = corr(XLOOPR,'Rows','pairwise');

fprintf('LOOPR X size: %s\n',mat2str(size(XLOOPR)));

%% ============================================================
% 3. IPIP-FMM
%
% Reconstruct the exact repeated samples used in the 100-run fit.
%% ============================================================

fmm_dir = fullfile(root_dir,'1. IPIP-FFM-data');

load(fullfile(fmm_dir,'Analysis','IPIPFMM_CC1.mat'), ...
    'IPIPFMM_CC1');

S = load(fullfile(fmm_dir,'IPIP_DDE_100_results6_comp.mat'));

assert(isfield(S,'repetition_seeds'), ...
    'repetition_seeds not found in IPIP-FMM results file.');

assert(isfield(S,'completed'), ...
    'completed not found in IPIP-FMM results file.');

assert(isfield(S,'sample_size'), ...
    'sample_size not found in IPIP-FMM results file.');

load(fullfile(fmm_dir,'B5_Red.mat'),'B5_Red');

item_names = [ ...
    "EXT"+string(1:10), ...
    "EST"+string(1:10), ...
    "AGR"+string(1:10), ...
    "CSN"+string(1:10), ...
    "OPN"+string(1:10)];

%% Complete-case population used for the comp fits

X_items = B5_Red{:,item_names};

complete_rows = all(ismember(X_items,1:5),2);

B5CM = B5_Red(complete_rows,:);

%% Same country collapsing used in the fitting workflow

country = categorical(B5CM.country);

country_names = categories(country);
country_counts = countcats(country);

large_countries = country_names(country_counts>=1000);

country_new = string(B5CM.country);
country_new(~ismember(country_new,string(large_countries))) = "OTHER";

B5CM.country = categorical(country_new);

fprintf('IPIP-FMM complete-case population: %d respondents\n',height(B5CM));

%% Final canonical item order used by averaged B1

row_order_fmm = IPIPFMM_CC1.item_map.OriginalRow(:)';

assert(numel(row_order_fmm)==50);
assert(isequal(sort(row_order_fmm),1:50));

%% Reconstruct repeated samples and average correlations

valid_iterations = find(S.completed);

fprintf('Completed IPIP-FMM repetitions: %d\n',numel(valid_iterations));

R_FMM_all = NaN(50,50,numel(valid_iterations));

for ii = 1:numel(valid_iterations)

    r = valid_iterations(ii);

    rng(S.repetition_seeds(r),'twister');

    Xr_table = rndIPIPFMM(B5CM,S.sample_size);

    % Extract raw 50 questionnaire responses
    Xr = Xr_table{:,item_names};

    % In case zero is used as missing
    Xr(Xr==0) = NaN;

    % Correlation first, in original EXT1,...,OPN10 order
    Rr = corr(Xr,'Rows','pairwise');

    % Apply exact final B1 item ordering to rows AND columns
    Rr = Rr(row_order_fmm,row_order_fmm);

    R_FMM_all(:,:,ii) = Rr;

    if(mod(ii,10)==0)
        fprintf('IPIP-FMM correlations completed: %d / %d\n', ...
            ii,numel(valid_iterations));
    end
end

RFMM = mean(R_FMM_all,3,'omitnan');

fprintf('IPIP-FMM averaged correlation size: %s\n',mat2str(size(RFMM)));

%% ============================================================
% 4. Sanity checks
%% ============================================================

assert(isequal(size(R100),[100 100]));
assert(isequal(size(RFMM),[50 50]));
assert(isequal(size(RLOOPR),[60 60]));

assert(all(abs(diag(R100)-1)<1e-10));
assert(all(abs(diag(RFMM)-1)<1e-10));
assert(all(abs(diag(RLOOPR)-1)<1e-10));

%% ============================================================
% 5. Common red-white-blue colormap
%% ============================================================

n_colors = 256;
n_half = floor(n_colors/2);

blue_white = [ ...
    linspace(0,1,n_half)', ...
    linspace(0,1,n_half)', ...
    ones(n_half,1)];

white_red = [ ...
    ones(n_colors-n_half,1), ...
    linspace(1,0,n_colors-n_half)', ...
    linspace(1,0,n_colors-n_half)'];

red_white_blue = [blue_white;white_red];

%% ============================================================
% 6. Plot
%% ============================================================

fig = figure('Color','w','Position',[80 100 1600 560]);

tl = tiledlayout(fig,1,3, ...
    'TileSpacing','compact', ...
    'Padding','compact');

ax1 = nexttile(tl);
plot_corr_panel( ...
    ax1, ...
    R100, ...
    string(IPIP100.domain_names), ...
    IPIP100.block_sizes, ...
    'IPIP-100');

ax2 = nexttile(tl);
plot_corr_panel( ...
    ax2, ...
    RFMM, ...
    string(IPIPFMM_CC1.domain_names), ...
    IPIPFMM_CC1.block_sizes, ...
    'IPIP-FFM');

ax3 = nexttile(tl);
plot_corr_panel( ...
    ax3, ...
    RLOOPR, ...
    string(LOOPR.domain_names), ...
    LOOPR.block_sizes, ...
    'LOOPR');

colormap(fig,red_white_blue);

cb = colorbar(ax3);
cb.Layout.Tile = 'east';
cb.Label.String = 'Correlation';
cb.FontSize = 12;

%% ============================================================
% 6. Plot
%% ============================================================

fig = figure('Color','w','Position',[80 100 1600 560]);

tl = tiledlayout(fig,1,3, ...
    'TileSpacing','compact', ...
    'Padding','compact');

ax1 = nexttile(tl);

plot_corr_panel( ...
    ax1, ...
    RFMM, ...
    string(IPIPFMM_CC1.domain_names), ...
    IPIPFMM_CC1.block_sizes, ...
    'IPIP-FFM');

ax2 = nexttile(tl);

plot_corr_panel( ...
    ax2, ...
    R100, ...
    string(IPIP100.domain_names), ...
    IPIP100.block_sizes, ...
    'IPIP-100');

ax3 = nexttile(tl);

plot_corr_panel( ...
    ax3, ...
    RLOOPR, ...
    string(LOOPR.domain_names), ...
    LOOPR.block_sizes, ...
    'LOOPR');

colormap(fig,red_white_blue);

cb = colorbar(ax3);
cb.Layout.Tile = 'east';
cb.Label.String = 'Correlation';
cb.FontSize = 12;
cb.Color = 'k';
cb.Label.Color = 'k';

%% ============================================================
% Local function
%% ============================================================

function plot_corr_panel(ax,R,domain_names,block_sizes,title_text)

imagesc(ax,R,[-1 1]);

axis(ax,'image');

set(ax, ...
    'YDir','reverse', ...
    'FontSize',12, ...
    'Box','on', ...
    'XColor','k', ...
    'YColor','k');

boundaries = cumsum(block_sizes);

starts = [1 boundaries(1:end-1)+1];

centers = (starts+boundaries)/2;

xticks(ax,centers);
yticks(ax,centers);

xticklabels(ax,domain_names);
yticklabels(ax,domain_names);

xlabel(ax,'Survey questions','Color','k');
ylabel(ax,'Survey questions','Color','k');

title(ax,title_text, ...
    'FontWeight','bold', ...
    'FontSize',16, ...
    'Color','k');

hold(ax,'on');

for j = 1:numel(boundaries)-1
    xline(ax,boundaries(j)+0.5,'k-','LineWidth',1.2);
    yline(ax,boundaries(j)+0.5,'k-','LineWidth',1.2);
end

hold(ax,'off');

end