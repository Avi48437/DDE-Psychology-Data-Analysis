function X_Red = rndIPIPFMM(B5CM, N)
%% Stratified sampling using country and elapsed-time groups

n_total = height(B5CM);

intro_group = discretize( ...
    B5CM.introelapse, ...
    quantile(B5CM.introelapse(~isnan(B5CM.introelapse)), [0 1/3 2/3 1]));
intro_group(isnan(intro_group)) = 4;

test_group = discretize( ...
    B5CM.testelapse, ...
    quantile(B5CM.testelapse(~isnan(B5CM.testelapse)), [0 1/3 2/3 1]));
test_group(isnan(test_group)) = 4;

end_group = discretize( ...
    B5CM.endelapse, ...
    quantile(B5CM.endelapse(~isnan(B5CM.endelapse)), [0 1/3 2/3 1]));
end_group(isnan(end_group)) = 4;

stratum = findgroups( ...
    B5CM.country, ...
    intro_group, ...
    test_group, ...
    end_group);

stratum_size = accumarray(stratum, 1);

sample_size = round(N * stratum_size / n_total);
sample_size = min(sample_size, stratum_size);

difference = N - sum(sample_size);

if difference > 0
    [~, order] = sort(stratum_size - sample_size, 'descend');
    sample_size(order(1:difference)) = ...
        sample_size(order(1:difference)) + 1;
elseif difference < 0
    [~, order] = sort(sample_size, 'descend');
    for k = 1:(-difference)
        sample_size(order(k)) = sample_size(order(k)) - 1;
    end
end

selected = false(n_total, 1);

for h = 1:length(stratum_size)

    idx = find(stratum == h);

    if sample_size(h) > 0
        idx = idx(randperm(length(idx), sample_size(h)));
        selected(idx) = true;
    end

end

X = B5CM(selected, :);

%% Keep only the 50 personality survey questions
survey_vars = [
    compose("EXT%d", 1:10), ...
    compose("EST%d", 1:10), ...
    compose("AGR%d", 1:10), ...
    compose("CSN%d", 1:10), ...
    compose("OPN%d", 1:10)
    ];

X_Red = X(:, survey_vars);

end