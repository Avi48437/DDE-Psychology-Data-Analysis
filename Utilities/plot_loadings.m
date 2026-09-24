function plot_loadings(B1_final, B2_final)

%% Separate intercept and loading columns

B1_intercept = B1_final(:, 1);
B1_loadings = B1_final(:, 2:end);

B2_intercept = B2_final(:, 1);
B2_loadings = B2_final(:, 2:end);

%% Check layer compatibility

if size(B1_loadings, 2) ~= size(B2_final, 1)
    error(['The number of B1 latent-factor columns must equal the ', ...
        'number of B2 rows.']);
end

%% Order B1 loading columns from strongest to weakest

B1_column_strength = vecnorm(B1_loadings, 2, 1);
[~, B1_column_order] = sort(B1_column_strength, 'descend');

B1_reordered = B1_loadings(:, B1_column_order);

%% Keep the B1 intercept as the leftmost column

B1_plot = [B1_intercept, B1_reordered];

%% Reorder all B2 rows using the B1 column order

B2_intercept = B2_intercept(B1_column_order, :);
B2_loadings = B2_loadings(B1_column_order, :);

%% Order only the non-intercept B2 columns

B2_column_strength = vecnorm(B2_loadings, 2, 1);
[~, B2_column_order] = sort(B2_column_strength, 'descend');

B2_reordered = B2_loadings(:, B2_column_order);

%% Keep the B2 intercept as the leftmost column

B2_plot = [B2_intercept, B2_reordered];

%% Use the same symmetric color scale

max_abs_loading = max(abs([B1_plot(:); B2_plot(:)]));

if max_abs_loading == 0
    max_abs_loading = 1;
end

color_limits = [-max_abs_loading, max_abs_loading];

%% Create red-white-blue colormap

number_of_colors = 256;
half_colors = floor(number_of_colors / 2);

blue_to_white = [
    linspace(0, 1, half_colors)', ...
    linspace(0, 1, half_colors)', ...
    ones(half_colors, 1)];

white_to_red = [
    ones(number_of_colors - half_colors, 1), ...
    linspace(1, 0, number_of_colors - half_colors)', ...
    linspace(1, 0, number_of_colors - half_colors)'];

red_white_blue = [blue_to_white; white_to_red];

%% Plot final B1 and B2

figure( ...
    'Position', [100 100 1400 650], ...
    'Color', 'w');

subplot(1, 2, 1)

imagesc(B1_plot)
set(gca, 'YDir', 'reverse')
axis tight
colorbar
caxis(color_limits)
colormap(red_white_blue)

title('Final B_1')
xlabel('Intercept and Layer 1 Factors')
ylabel('Survey Items')

xticks(1:size(B1_plot, 2))
xticklabels(["Intercept", "A1_" + string(B1_column_order)])

yticks(1:size(B1_plot, 1))
yticklabels(1:size(B1_plot, 1))

subplot(1, 2, 2)

imagesc(B2_plot)
set(gca, 'YDir', 'reverse')
axis tight
colorbar
caxis(color_limits)
colormap(red_white_blue)

title('Final B_2')
xlabel('Intercept and Layer 2 Factors')
ylabel('Reordered Layer 1 Factors')

xticks(1:size(B2_plot, 2))
xticklabels(["Intercept", "A2_" + string(B2_column_order)])

yticks(1:size(B2_plot, 1))
yticklabels(B1_column_order)

sgtitle('Final DDE Loading Matrices')

end