function plot_loadings2(DDE,block_size)

B1 = DDE.B1;
B2 = DDE.B2;

%% ============================================================
% Basic checks
%% ============================================================

if size(B1,2)-1 ~= size(B2,1)
    error('Number of B1 factors must equal number of B2 rows.');
end

if mod(size(B1,1),block_size) ~= 0
    error('Number of B1 rows must be divisible by block_size.');
end

n_domains = size(B1,1)/block_size;

if isfield(DDE,'domain_names')
    domain_names = string(DDE.domain_names);
else
    domain_names = ["EXT","EST","AGR","CSN","OPN"];
end

if numel(domain_names) ~= n_domains
    error('Number of domain names does not match the number of blocks.');
end

%% ============================================================
% No reordering here
%% ============================================================

B1_plot = B1;
B2_plot = B2;

n_A1 = size(B1_plot,2)-1;
n_A2 = size(B2_plot,2)-1;

%% ============================================================
% Common symmetric color scale
%% ============================================================

max_abs_loading = max(abs([B1_plot(:);B2_plot(:)]));

if max_abs_loading == 0
    max_abs_loading = 1;
end

color_limits = [-max_abs_loading max_abs_loading];

%% ============================================================
% Red-white-blue colormap
%% ============================================================

number_of_colors = 256;
half_colors = floor(number_of_colors/2);

blue_to_white = [linspace(0,1,half_colors)',linspace(0,1,half_colors)',ones(half_colors,1)];
white_to_red = [ones(number_of_colors-half_colors,1),linspace(1,0,number_of_colors-half_colors)',linspace(1,0,number_of_colors-half_colors)'];
red_white_blue = [blue_to_white;white_to_red];

%% ============================================================
% Domain positions
%% ============================================================

domain_centers = block_size/2 + 0.5 + (0:n_domains-1)*block_size;
domain_boundaries = block_size + 0.5:block_size:size(B1_plot,1);

%% ============================================================
% Figure
%% ============================================================

figure('Position',[100 100 1500 700],'Color','w');

%% ============================================================
% B1
%% ============================================================

subplot(1,2,1)

imagesc(B1_plot)

set(gca,'YDir','reverse')
axis tight
caxis(color_limits)
colormap(red_white_blue)

cb1 = colorbar;
cb1.Color = 'k';
cb1.FontSize = 12;
cb1.LineWidth = 1;

title('B_1','FontSize',18,'FontWeight','bold','Color','k')
xlabel('Intercept and Layer 1 Factors','FontSize',16,'FontWeight','normal','Color','k')
ylabel('Big Five Domain','FontSize',16,'FontWeight','normal','Color','k')

xticks(1:size(B1_plot,2))
xticklabels(["Intercept","A1_" + string(1:n_A1)])
xtickangle(35)

yticks(domain_centers)
yticklabels(domain_names)

ax1 = gca;
ax1.FontSize = 14;
ax1.FontWeight = 'normal';
ax1.XColor = 'k';
ax1.YColor = 'k';
ax1.LineWidth = 1;
ax1.TickLength = [0.005 0.005];
ax1.YAxis.TickLabelRotation = 90;

hold on
for y = domain_boundaries
    yline(y,'k-','LineWidth',1.2);
end
hold off

%% ============================================================
% B2
%% ============================================================

subplot(1,2,2)

imagesc(B2_plot)

set(gca,'YDir','reverse')
axis tight
caxis(color_limits)
colormap(red_white_blue)

cb2 = colorbar;
cb2.Color = 'k';
cb2.FontSize = 12;
cb2.LineWidth = 1;

title('B_2','FontSize',18,'FontWeight','bold','Color','k')
xlabel('Intercept and Layer 2 Factors','FontSize',16,'FontWeight','normal','Color','k')
ylabel('Layer 1 Factors','FontSize',16,'FontWeight','normal','Color','k')

xticks(1:size(B2_plot,2))
xticklabels(["Intercept","A2_" + string(1:n_A2)])
xtickangle(0)

yticks(1:n_A1)
yticklabels("A1_" + string(1:n_A1))

ax2 = gca;
ax2.FontSize = 14;
ax2.FontWeight = 'normal';
ax2.XColor = 'k';
ax2.YColor = 'k';
ax2.LineWidth = 1;
ax2.TickLength = [0.005 0.005];

%% ============================================================
% Overall title
%% ============================================================

if isfield(DDE,'dataset_name')
    sgtitle(string(DDE.dataset_name) + " DDE Loading Matrices",'FontSize',22,'FontWeight','bold','Color','k')
else
    sgtitle('DDE Loading Matrices','FontSize',22,'FontWeight','bold','Color','k')
end

end