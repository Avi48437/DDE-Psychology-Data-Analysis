%% ============================================================
% Paths
%% ============================================================

root_dir = '/Users/avinandanroy/Desktop/Research/DDE/DDE_Fitting_Pshycology2';

addpath(fullfile(root_dir,'Scripts'));
addpath(fullfile(root_dir,'Algorithms'));
addpath(fullfile(root_dir,'Utilities'));

ipipfmm_file = fullfile( ...
    root_dir, ...
    '1. IPIP-FFM-data', ...
    'Analysis', ...
    'PostProcessed_Data', ...
    'IPIPFMM_CC1.mat');

ipip100_file = fullfile( ...
    root_dir, ...
    '2. IPIP 100', ...
    'Analysis', ...
    'PostProcessed_Data', ...
    'IPIP100_PostProcessed.mat');

ipip98_file = fullfile( ...
    root_dir, ...
    '2. IPIP 100', ...
    'Analysis', ...
    'PostProcessed_Data', ...
    'IPIP98_PostProcessed.mat');

loopr_file = fullfile( ...
    root_dir, ...
    'LOPR', ...
    'Analysis', ...
    'PostProcessed_Data', ...
    'LOOPR_PostProcessed.mat');

%% ============================================================
% Load Fits
%% ============================================================
IPIPFMM_CC1 = load(ipipfmm_file).IPIPFMM_CC1;
IPIP100 = load(ipip100_file).IPIP100;
IPIP98 = load(ipip98_file).IPIP98;
LOOPR = load(loopr_file).LOOPR;
%% ============================================================
% Loaded-model summary
%% ============================================================

fprintf('\n============================================================\n');
fprintf('POST-PROCESSED DDE FITS LOADED\n');
fprintf('============================================================\n');
fprintf('IPIP-FMM C1 : J = %d, K1 = %d, K2 = %d\n',size(IPIPFMM_CC1.B1,1),IPIPFMM_CC1.K1,IPIPFMM_CC1.K2);
fprintf('IPIP-100    : J = %d, K1 = %d, K2 = %d\n',size(IPIP100.B1,1),IPIP100.K1,IPIP100.K2);
fprintf('IPIP-98    : J = %d, K1 = %d, K2 = %d\n',size(IPIP100.B1,1),IPIP100.K1,IPIP100.K2);
fprintf('LOOPR BFI-2 : J = %d, K1 = %d, K2 = %d\n',size(LOOPR.B1,1),LOOPR.K1,LOOPR.K2);
fprintf('============================================================\n');

%% ============================================================
% Threshold IPIP-FMM C1
%
% Keep Layer-1 factors whose relative column L1 norm is at
% least the specified threshold.
%% ============================================================

threshold = 0.40;
IPIPFMM = threshold_model(IPIPFMM_CC1,threshold);
IPIPFMM.dataset_name = "IPIP-FMM";

%% ============================================================
% Plot loading matrices
%
% true:
%   B1 | B2 | B1 clustered by Layer 2 | B2 clustered
%% ============================================================
set(groot,'defaultFigureWindowStyle','docked');
plot_loadings_mult( ...
    {'IPIPFMM_CC1','IPIP98','LOOPR'});
plot_loadings_mult( ...
    {'IPIPFMM','IPIP98','LOOPR'});

% ============================================================
% Hierarchical anchor maps
% ============================================================

anchor_map(IPIPFMM_CC1,3,'network')
anchor_map(IPIPFMM,3,'network');
anchor_map(IPIP100,3,'network');
anchor_map(IPIP98,3,'network');
anchor_map(LOOPR,3,'network');