function plot_loadings_comparison(B1_FMM,B2_FMM,B1_IPIP100,B2_IPIP100)

%% ============================================================
% Active Layer-1 factors
% ============================================================

active_FMM = find(any(B1_FMM(:,2:end) ~= 0,1));
active_100 = find(any(B1_IPIP100(:,2:end) ~= 0,1));

%% ============================================================
% Arrange Layer-1 factors by personality domain
%
% IPIP-FMM:
% EXT : (1,14)
% EST : (2,12)
% AGR : (4,15)
% CSN : (6,3)
% OPN : (11,7,5)
%
% Any remaining active factor is appended afterward.
% ============================================================

desired_FMM = [1 14 2 12 4 15 6 3 11 7 5];

order_FMM = desired_FMM(ismember(desired_FMM,active_FMM));
remaining_FMM = active_FMM(~ismember(active_FMM,order_FMM));
order_FMM = [order_FMM remaining_FMM];

%% ============================================================
% IPIP-100 Layer-1 ordering
%
% EXT : (1,2)
% EST : (8,18)
% AGR : (13)
% CSN : (5,32)
% OPN : (10)
% ============================================================

desired_100 = [1 2 8 18 13 5 32 10];

order_100 = desired_100(ismember(desired_100,active_100));
remaining_100 = active_100(~ismember(active_100,order_100));
order_100 = [order_100 remaining_100];

%% ============================================================
% Reorder IPIP-FMM survey-item rows WITHIN each domain
%
% Each domain contains 10 rows.
% Sorting factor:
% EXT -> A1_1
% EST -> A1_2
% AGR -> A1_4
% CSN -> A1_6
% OPN -> A1_11
% ============================================================

domain_blocks_FMM = {1:10,11:20,21:30,31:40,41:50};
anchor_factors_FMM = [1 2 4 6 11];

row_order_FMM = zeros(1,50);
pos = 1;

for d = 1:5

    rows = domain_blocks_FMM{d};
    anchor = anchor_factors_FMM(d);

    loading = B1_FMM(rows,anchor+1);
    [~,idx] = sort(loading,'descend');

    row_order_FMM(pos:pos+numel(rows)-1) = rows(idx);
    pos = pos + numel(rows);

end

%% ============================================================
% Reorder IPIP-100 survey-item rows WITHIN each domain
%
% B1_IPIP100 is already grouped:
% 1:20   EXT
% 21:40  EST
% 41:60  AGR
% 61:80  CSN
% 81:100 OPN
%
% Sorting factor:
% EXT -> A1_1
% EST -> A1_8
% AGR -> A1_13
% CSN -> A1_5
% OPN -> A1_10
% ============================================================

domain_blocks_100 = {1:20,21:40,41:60,61:80,81:100};
anchor_factors_100 = [1 8 13 5 10];

row_order_100 = zeros(1,100);
pos = 1;

for d = 1:5

    rows = domain_blocks_100{d};
    anchor = anchor_factors_100(d);

    loading = B1_IPIP100(rows,anchor+1);
    [~,idx] = sort(loading,'descend');

    row_order_100(pos:pos+numel(rows)-1) = rows(idx);
    pos = pos + numel(rows);

end

%% ============================================================
% Construct B1 matrices
% Keep intercept + active ordered Layer-1 factors
% ============================================================

B1_plot_FMM = [ ...
    B1_FMM(row_order_FMM,1), ...
    B1_FMM(row_order_FMM,order_FMM+1)];

B1_plot_100 = [ ...
    B1_IPIP100(row_order_100,1), ...
    B1_IPIP100(row_order_100,order_100+1)];

%% ============================================================
% Construct B2 matrices
%
% B2 rows correspond to Layer-1 factors.
% Therefore use exactly the same Layer-1 ordering as B1 columns.
% ============================================================

B2_rows_FMM = B2_FMM(order_FMM,:);
B2_rows_100 = B2_IPIP100(order_100,:);

%% Keep only active Layer-2 factors, but retain intercept

active2_FMM = find(any(B2_rows_FMM(:,2:end) ~= 0,1));
active2_100 = find(any(B2_rows_100(:,2:end) ~= 0,1));

B2_plot_FMM = [ ...
    B2_rows_FMM(:,1), ...
    B2_rows_FMM(:,active2_FMM+1)];

B2_plot_100 = [ ...
    B2_rows_100(:,1), ...
    B2_rows_100(:,active2_100+1)];

%% ============================================================
% Common symmetric color scale across all four panels
% ============================================================

max_abs_loading = max(abs([ ...
    B1_plot_FMM(:); ...
    B2_plot_FMM(:); ...
    B1_plot_100(:); ...
    B2_plot_100(:)]));

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
% Plot
% ============================================================

figure('Position',[50 50 1600 1050],'Color','w');

tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

%% ============================================================
% IPIP-FMM B1
% ============================================================

nexttile

imagesc(B1_plot_FMM)
set(gca,'YDir','reverse')
axis tight
caxis(color_limits)
colormap(red_white_blue)
colorbar

title('IPIP-FMM C1: B_1  (12,1)')
xlabel('Intercept and Layer 1 Factors')
ylabel('Personality Trait')

xticks(1:size(B1_plot_FMM,2))
xticklabels(["Intercept","A1_" + string(order_FMM)])
xtickangle(45)

% One label per 10-item personality block
yticks([5.5 15.5 25.5 35.5 45.5])
yticklabels({'EXT','EST','AGR','CSN','OPN'})

% Domain separators
yline(10.5,'k-','LineWidth',1.2)
yline(20.5,'k-','LineWidth',1.2)
yline(30.5,'k-','LineWidth',1.2)
yline(40.5,'k-','LineWidth',1.2)

%% ============================================================
% IPIP-FMM B2
% ============================================================

nexttile

imagesc(B2_plot_FMM)
set(gca,'YDir','reverse')
axis tight
caxis(color_limits)
colormap(red_white_blue)
colorbar

title('IPIP-FMM C1: B_2  (12,1)')
xlabel('Intercept and Layer 2 Factors')
ylabel('Layer 1 Factors')

xticks(1:size(B2_plot_FMM,2))
xticklabels(["Intercept","A2_" + string(active2_FMM)])

yticks(1:numel(order_FMM))
yticklabels("A1_" + string(order_FMM))

%% ============================================================
% IPIP-100 B1
% ============================================================

nexttile

imagesc(B1_plot_100)
set(gca,'YDir','reverse')
axis tight
caxis(color_limits)
colormap(red_white_blue)
colorbar

title('IPIP-100: B_1  (8,1)')
xlabel('Intercept and Layer 1 Factors')
ylabel('Personality Trait')

xticks(1:size(B1_plot_100,2))
xticklabels(["Intercept","A1_" + string(order_100)])
xtickangle(45)

% One label per 20-item personality block
yticks([10.5 30.5 50.5 70.5 90.5])
yticklabels({'EXT','EST','AGR','CSN','OPN'})

% Domain separators
yline(20.5,'k-','LineWidth',1.2)
yline(40.5,'k-','LineWidth',1.2)
yline(60.5,'k-','LineWidth',1.2)
yline(80.5,'k-','LineWidth',1.2)

%% ============================================================
% IPIP-100 B2
% ============================================================

nexttile

imagesc(B2_plot_100)
set(gca,'YDir','reverse')
axis tight
caxis(color_limits)
colormap(red_white_blue)
colorbar

title('IPIP-100: B_2  (8,1)')
xlabel('Intercept and Layer 2 Factors')
ylabel('Layer 1 Factors')

xticks(1:size(B2_plot_100,2))
xticklabels(["Intercept","A2_" + string(active2_100)])

yticks(1:numel(order_100))
yticklabels("A1_" + string(order_100))

%% ============================================================
% Overall title
% ============================================================

title(tl,'DDE Loading Matrices: IPIP-FMM C1 vs IPIP-100')

end