clear
clc
close all

%% ============================================================
% Load best LOOPR DDE model selected by pMSE tuning
%% ============================================================

script_dir = fileparts(mfilename('fullpath'));
loopr_dir = fileparts(fileparts(script_dir));

best_model_file = fullfile(loopr_dir,'Fitting','Results','LOOPR_best_DDE_tuning.mat');


if ~isfile(best_model_file)
    error('Best-model file not found: %s',best_model_file);
end

best_fit = load(best_model_file);

B1 = best_fit.best_B1_CSP;
B2 = best_fit.best_B2_CSP;

num_act1 = best_fit.best_num_act1;
num_act2 = best_fit.best_num_act2;

best_settings = best_fit.best_settings;
best_pMSE = best_fit.best_avg_pMSE;

fprintf('\n============================================================\n');
fprintf('Best LOOPR DDE model\n');
fprintf('============================================================\n');
fprintf('tau  = %.4f\n',best_settings.tau);
fprintf('t1   = %.4f\n',best_settings.t1);
fprintf('t2   = %.4f\n',best_settings.t2);
fprintf('temp = %.4f\n',best_settings.temp);
fprintf('Active B1 factors = %d\n',num_act1);
fprintf('Active B2 factors = %d\n',num_act2);
fprintf('Mean pMSE = %.8f\n',best_pMSE);
fprintf('============================================================\n');

%% ============================================================
% Group BFI-2 survey items by domain
%
% EXT : 1,6,11,...,56
% AGR : 2,7,12,...,57
% CSN : 3,8,13,...,58
% NEM : 4,9,14,...,59
% OPN : 5,10,15,...,60
%% ============================================================

survey_order = [ ...
    1:5:60, ...
    2:5:60, ...
    3:5:60, ...
    4:5:60, ...
    5:5:60];

domain_names = ["EXT","AGR","CSN","NEM","OPN"];

B1 = B1(survey_order,:);

%% ============================================================
% Separate intercept and loading matrices
%% ============================================================

B1_intercept = B1(:,1);
B1_loadings = B1(:,2:end);

B2_intercept = B2(:,1);
B2_loadings = B2(:,2:end);

if size(B1_loadings,2) ~= size(B2,1)
    error('The number of B1 factors must equal the number of B2 rows.');
end

%% ============================================================
% Identify active Layer-1 factors
%% ============================================================

B1_strength = vecnorm(B1_loadings,2,1);

[~,strength_order] = sort(B1_strength,'descend');

if num_act1 > length(strength_order)
    error('num_act1 exceeds the number of available B1 factors.');
end

B1_active = strength_order(1:num_act1);

%% ============================================================
% Compute domain strength of each Layer-1 factor
%
% For each factor, calculate the L2 norm of its 12 loadings
% within EXT, AGR, CSN, NEM and OPN.
%% ============================================================

domain_strength = zeros(5,size(B1_loadings,2));

for d = 1:5
    rows = (d-1)*12 + (1:12);
    domain_strength(d,:) = vecnorm(B1_loadings(rows,:),2,1);
end

%% ============================================================
% Assign each active B1 factor to its strongest domain
%% ============================================================

[~,anchor_domain] = max(domain_strength,[],1);

%% ============================================================
% Permute active B1 columns by domain
%
% Domain order:
% EXT -> AGR -> CSN -> NEM -> OPN
%
% Within each domain:
% strongest domain factor comes first.
%% ============================================================

B1_column_order = [];
domain_counts = zeros(1,5);

for d = 1:5

    idx = B1_active(anchor_domain(B1_active) == d);

    if ~isempty(idx)
        [~,ord] = sort(domain_strength(d,idx),'descend');
        idx = idx(ord);
    end

    domain_counts(d) = numel(idx);

    B1_column_order = [B1_column_order,idx];

end

if length(B1_column_order) ~= num_act1
    error('Not all active B1 factors were assigned to a domain.');
end

%% ============================================================
% Construct B1 with final COLUMN permutation
%% ============================================================

B1_plot = [B1_intercept,B1_loadings(:,B1_column_order)];

%% ============================================================
% Apply EXACTLY the same Layer-1 factor permutation to B2 rows
%% ============================================================

B2_intercept = B2_intercept(B1_column_order,:);
B2_loadings = B2_loadings(B1_column_order,:);

%% ============================================================
% Permute survey rows WITHIN EACH DOMAIN
%
% For each domain:
%
% 1. Find the FIRST displayed B1 factor belonging to that domain.
% 2. Take its 12 loading values.
% 3. Sort those values from HIGHEST to LOWEST.
% 4. Move the ENTIRE B1 rows using that same permutation.
%
% Therefore:
%
% strongest positive/red -> TOP
% near zero/white        -> MIDDLE
% strongest negative/blue-> BOTTOM
%% ============================================================

final_survey_order = survey_order;

for d = 1:5

    rows = (d-1)*12 + (1:12);

    if domain_counts(d) == 0
        continue
    end

    % First displayed factor column for this domain.
    % Column 1 of B1_plot is the intercept.
    reference_column = 2 + sum(domain_counts(1:d-1));

    % Actual 12 loading values
    reference_values = B1_plot(rows,reference_column);

    % Highest -> lowest
    [~,row_perm] = sort(reference_values,'descend');

    % Save entire current block
    current_block = B1_plot(rows,:);

    % Permute entire rows
    B1_plot(rows,:) = current_block(row_perm,:);

    % Keep track of original BFI item identities
    current_order = final_survey_order(rows);
    final_survey_order(rows) = current_order(row_perm);

end

%% ============================================================
% Verify domain-wise row permutation
%
% First factor of every domain must decrease monotonically
% from top to bottom.
%% ============================================================

for d = 1:5

    if domain_counts(d) == 0
        continue
    end

    rows = (d-1)*12 + (1:12);

    reference_column = 2 + sum(domain_counts(1:d-1));

    vals = B1_plot(rows,reference_column);

    if any(diff(vals) > 1e-12)
        error('Row permutation failed for domain %s.',domain_names(d));
    end

end

fprintf('\nDomain-wise B1 row ordering verified successfully.\n');

%% ============================================================
% Identify active Layer-2 factors
%% ============================================================

B2_strength = vecnorm(B2_loadings,2,1);

[~,B2_strength_order] = sort(B2_strength,'descend');

if num_act2 > length(B2_strength_order)
    error('num_act2 exceeds the number of available B2 factors.');
end

B2_column_order = B2_strength_order(1:num_act2);

%% ============================================================
% Construct B2 plotting matrix
%% ============================================================

B2_plot = [B2_intercept,B2_loadings(:,B2_column_order)];

%% ============================================================
% Axis labels
%% ============================================================

B1_labels = "A1_" + string(B1_column_order);
B2_labels = "A2_" + string(B2_column_order);

%% ============================================================
% Domain boundaries among B1 factors
%% ============================================================

domain_boundaries = cumsum(domain_counts);

%% ============================================================
% Common symmetric color scale
%% ============================================================

max_abs_loading = max(abs([B1_plot(:);B2_plot(:)]));

if max_abs_loading == 0
    max_abs_loading = 1;
end

color_limits = [-max_abs_loading,max_abs_loading];

%% ============================================================
% Red-white-blue colormap
%% ============================================================

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

red_white_blue = [blue_to_white;white_to_red];

%% ============================================================
% Plot
%% ============================================================

figure('Position',[100 100 1500 700],'Color','w');

% ============================================================
% B1
% ============================================================

subplot(1,2,1)

imagesc(B1_plot)

set(gca,'YDir','reverse')
axis tight

colorbar
caxis(color_limits)
colormap(red_white_blue)

hold on

% Horizontal Big Five domain boundaries
yline(12.5,'k-','LineWidth',1.2)
yline(24.5,'k-','LineWidth',1.2)
yline(36.5,'k-','LineWidth',1.2)
yline(48.5,'k-','LineWidth',1.2)

% Vertical Layer-1 domain boundaries
for d = 1:4

    if domain_boundaries(d) > 0 && ...
       domain_boundaries(d) < length(B1_column_order)

        xline(domain_boundaries(d)+1.5,'k-','LineWidth',1.2);

    end

end

hold off

title('LOOPR: B_1')

xlabel('Intercept and Layer 1 Factors')
ylabel('Personality Trait')

xticks(1:size(B1_plot,2))
xticklabels(["Intercept",B1_labels])
xtickangle(45)

yticks([6.5 18.5 30.5 42.5 54.5])
yticklabels(domain_names)

% ============================================================
% B2
% ============================================================

subplot(1,2,2)

imagesc(B2_plot)

set(gca,'YDir','reverse')
axis tight

colorbar
caxis(color_limits)
colormap(red_white_blue)

hold on

% Apply same Layer-1 domain grouping to B2 rows
for d = 1:4

    if domain_boundaries(d) > 0 && ...
       domain_boundaries(d) < length(B1_column_order)

        yline(domain_boundaries(d)+0.5,'k-','LineWidth',1.2);

    end

end

hold off

title('LOOPR: B_2')

xlabel('Intercept and Layer 2 Factors')
ylabel('Layer 1 Factors')

xticks(1:size(B2_plot,2))
xticklabels(["Intercept",B2_labels])
xtickangle(45)

yticks(1:length(B1_column_order))
yticklabels(B1_labels)

%% ============================================================
% Overall title
%% ============================================================

sgtitle(sprintf( ...
    'LOOPR DDE | tau=%.2f, t1=%.3f, t2=%.3f, temp=%.2f, pMSE=%.6f', ...
    best_settings.tau, ...
    best_settings.t1, ...
    best_settings.t2, ...
    best_settings.temp, ...
    best_pMSE));