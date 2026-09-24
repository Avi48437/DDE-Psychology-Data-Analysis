%% Group IPIP-100 rows by Big-Five domain
% Order chosen to match current plots:
% EXT | EST/NEU | AGR | CSN | OPN

EXT_idx = [3 10 14 18 23 29 33 39 43 49 53 59 63 69 73 79 83 89 93 99];

EST_idx = [11 12 17 19 27 30 37 40 47 50 57 60 67 70 77 80 87 90 97 100];

AGR_idx = [2 6 9 13 22 26 32 36 42 46 52 56 62 66 72 76 82 86 92 96];

CSN_idx = [5 8 15 20 25 28 35 38 45 48 55 58 65 68 75 78 85 88 95 98];

OPN_idx = [1 4 7 16 21 24 31 34 41 44 51 54 61 64 71 74 81 84 91 94];

row_order = [EXT_idx EST_idx AGR_idx CSN_idx OPN_idx];

%% Sanity check
fprintf('Number of reordered rows: %d\n',numel(row_order));
fprintf('All rows represented exactly once: %d\n', ...
    isequal(sort(row_order),1:100));

%% Reorder only B1 rows
B1_IPIP100_grouped = results.B1_CSP(row_order,:);
B2_IPIP100 = results.B2_CSP;

%% Plot
plot_loadings(B1_IPIP100_grouped,B2_IPIP100);
