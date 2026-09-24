clear
clc
close all

%% ============================================================
% IPIP-100 / IPIP-98 Predictive Analysis Data
%
% Conventional Big-Five scoring source:
%
% Official IPIP NEO-domain scoring key:
% https://www.ipip.ori.org/newNEODomainsKey.html
%
% The original BIGFIVE responses are coded on a 1,...,5 scale.
% Negatively keyed items are reverse coded as
%
%                       x_rev = 6 - x.
%
% After reverse coding, the conventional score for each domain
% is the respondent-wise average of all items belonging to that
% domain:
%
%               EXT | EST | AGR | CSN | OPN
%
% Missing responses, if present, are omitted from the average.
%
% IPIP-100:
%   Uses all 100 BIGFIVE items, with 20 items per domain.
%
% IPIP-98:
%   Uses the same questionnaire but removes:
%
%       BIGFIVE_44 : conservative political-candidate item
%       BIGFIVE_51 : liberal political-candidate item
%
%   Both removed items belong to OPN. Therefore IPIP-98 uses
%   20 EXT, 20 EST, 20 AGR, 20 CSN, and 18 OPN items.
%
% Saved data:
%
% rawX100 : original 100 questionnaire responses
% rawX98  : original responses excluding BIGFIVE_44 and BIGFIVE_51
%
% X100    : five reverse-coded IPIP-100 domain averages
% X98     : five reverse-coded IPIP-98 domain averages
%
% A100    : canonical IPIP-100 DDE coordinates [A1,A2]
% A98     : canonical IPIP-98 DDE coordinates [A1,A2]
%
% Y       : common downstream outcomes
%% ============================================================

%% Paths

root_dir = '/Users/avinandanroy/Desktop/Research/DDE/DDE_Fitting_Pshycology2';

ipip_dir = fullfile(root_dir,'2. IPIP 100');
data_dir = fullfile(ipip_dir,'Data');
analysis_dir = fullfile(ipip_dir,'Analysis','PostProcessed_Data');
predictive_dir = fullfile(ipip_dir,'Analysis','Predictive_Data');

if ~exist(predictive_dir,'dir')
    mkdir(predictive_dir);
end

%% ============================================================
% 1. Load post-processed DDE representations
%% ============================================================

load(fullfile(analysis_dir,'IPIP100_PostProcessed.mat'),'IPIP100');
load(fullfile(analysis_dir,'IPIP98_PostProcessed.mat'),'IPIP98');

%% ============================================================
% 2. Load original enriched questionnaire data
%% ============================================================

B5_wo = readtable(fullfile(data_dir,'B5_wo.csv'));

N = height(B5_wo);
available_names = string(B5_wo.Properties.VariableNames);

fprintf('\nB5_wo size: %d x %d\n',height(B5_wo),width(B5_wo));

%% ============================================================
% 3. Check respondent dimensions
%% ============================================================

assert(size(IPIP100.A1,1) == N,'IPIP100 A1 rows do not match B5_wo.');
assert(size(IPIP100.A2,1) == N,'IPIP100 A2 rows do not match B5_wo.');

assert(size(IPIP98.A1,1) == N,'IPIP98 A1 rows do not match B5_wo.');
assert(size(IPIP98.A2,1) == N,'IPIP98 A2 rows do not match B5_wo.');

%% ============================================================
% 4. Extract original BIGFIVE_1,...,BIGFIVE_100 responses
%% ============================================================

big5_names = "BIGFIVE_" + string(1:100);

assert(all(ismember(big5_names,available_names)), ...
    'Some BIGFIVE_1,...,BIGFIVE_100 variables are missing.');

rawX100_matrix = B5_wo{:,big5_names};

assert(size(rawX100_matrix,2) == 100);

rawX100 = array2table(rawX100_matrix,'VariableNames',cellstr(big5_names));

%% ============================================================
% 5. Construct raw IPIP-98 questionnaire matrix
%
% Remove original questions 44 and 51.
%% ============================================================

q98 = setdiff(1:100,[44 51],'stable');

rawX98_matrix = rawX100_matrix(:,q98);
rawX98_names = "BIGFIVE_" + string(q98);

rawX98 = array2table(rawX98_matrix,'VariableNames',cellstr(rawX98_names));

assert(width(rawX98) == 98);
assert(~ismember("BIGFIVE_44",string(rawX98.Properties.VariableNames)));
assert(~ismember("BIGFIVE_51",string(rawX98.Properties.VariableNames)));

%% ============================================================
% 6. Big-Five domain membership
%% ============================================================

EXT_idx = [3 10 14 18 23 29 33 39 43 49 53 59 63 69 73 79 83 89 93 99];

EST_idx = [11 12 17 19 27 30 37 40 47 50 57 60 67 70 77 80 87 90 97 100];

AGR_idx = [2 6 9 13 22 26 32 36 42 46 52 56 62 66 72 76 82 86 92 96];

CSN_idx = [5 8 15 20 25 28 35 38 45 48 55 58 65 68 75 78 85 88 95 98];

OPN_idx = [1 4 7 16 21 24 31 34 41 44 51 54 61 64 71 74 81 84 91 94];

%% ============================================================
% 7. Negatively keyed items
%
% Reverse coding:
%
%                     x_rev = 6 - x
%% ============================================================

reverse_q = [ ...
    14 18 29 39 49 59 69 79 89 99, ...
    12 17 30 40 50 60 70 80 90 100, ...
    2 9 22 32 42 52 62 72 82 92, ...
    8 20 28 38 48 58 68 78 88 98, ...
    4 7 24 34 44 54 64 74 84 94];

%% ============================================================
% 8. Reverse-code the 100-item questionnaire
%% ============================================================

X_scored = rawX100_matrix;

X_scored(:,reverse_q) = 6-X_scored(:,reverse_q);

%% ============================================================
% 9. Construct conventional IPIP-100 Big-Five scores
%
% Each score is the mean of the 20 reverse-coded items
% belonging to that domain.
%% ============================================================

EXT100 = mean(X_scored(:,EXT_idx),2,'omitnan');
EST100 = mean(X_scored(:,EST_idx),2,'omitnan');
AGR100 = mean(X_scored(:,AGR_idx),2,'omitnan');
CSN100 = mean(X_scored(:,CSN_idx),2,'omitnan');
OPN100 = mean(X_scored(:,OPN_idx),2,'omitnan');

X100 = table(EXT100,EST100,AGR100,CSN100,OPN100, ...
    'VariableNames',{'EXT','EST','AGR','CSN','OPN'});

%% ============================================================
% 10. Construct conventional IPIP-98 Big-Five scores
%
% q44 and q51 are omitted.
%
% EXT, EST, AGR, and CSN therefore remain 20-item averages.
% OPN becomes an 18-item average.
%% ============================================================

OPN_idx_98 = setdiff(OPN_idx,[44 51],'stable');

assert(numel(EXT_idx) == 20);
assert(numel(EST_idx) == 20);
assert(numel(AGR_idx) == 20);
assert(numel(CSN_idx) == 20);
assert(numel(OPN_idx_98) == 18);

EXT98 = mean(X_scored(:,EXT_idx),2,'omitnan');
EST98 = mean(X_scored(:,EST_idx),2,'omitnan');
AGR98 = mean(X_scored(:,AGR_idx),2,'omitnan');
CSN98 = mean(X_scored(:,CSN_idx),2,'omitnan');
OPN98 = mean(X_scored(:,OPN_idx_98),2,'omitnan');

X98 = table(EXT98,EST98,AGR98,CSN98,OPN98, ...
    'VariableNames',{'EXT','EST','AGR','CSN','OPN'});

%% ============================================================
% 11. Construct IPIP-100 DDE representation A100
%% ============================================================

K1_100 = size(IPIP100.A1,2);
K2_100 = size(IPIP100.A2,2);

A100_matrix = [IPIP100.A1 IPIP100.A2];

A100_names = [ ...
    "A1_" + string(1:K1_100), ...
    "A2_" + string(1:K2_100)];

A100 = array2table(A100_matrix,'VariableNames',cellstr(A100_names));

%% ============================================================
% 12. Construct IPIP-98 DDE representation A98
%% ============================================================

K1_98 = size(IPIP98.A1,2);
K2_98 = size(IPIP98.A2,2);

A98_matrix = [IPIP98.A1 IPIP98.A2];

A98_names = [ ...
    "A1_" + string(1:K1_98), ...
    "A2_" + string(1:K2_98)];

A98 = array2table(A98_matrix,'VariableNames',cellstr(A98_names));

%% ============================================================
% 13. Construct common outcome matrix Y
%% ============================================================

outcome_names = [ ...
    "educ", ...
    "faminc_new", ...
    "marstat", ...
    "pew_prayer", ...
    "votereg", ...
    "pid3", ...
    "pid7", ...
    "ideo5"];

assert(all(ismember(outcome_names,available_names)), ...
    'One or more requested outcomes are missing from B5_wo.');

Y = B5_wo(:,outcome_names);

%% ============================================================
% 14. Final alignment checks
%% ============================================================

assert(height(rawX100) == N);
assert(height(rawX98) == N);

assert(height(X100) == N);
assert(height(X98) == N);

assert(height(A100) == N);
assert(height(A98) == N);

assert(height(Y) == N);

assert(width(rawX100) == 100);
assert(width(rawX98) == 98);

assert(width(X100) == 5);
assert(width(X98) == 5);

assert(width(A100) == K1_100+K2_100);
assert(width(A98) == K1_98+K2_98);

%% ============================================================
% 15. Summary
%% ============================================================

fprintf('\n============================================================\n');
fprintf('IPIP PREDICTIVE ANALYSIS DATA\n');
fprintf('============================================================\n');

fprintf('rawX100 : %d x %d\n',height(rawX100),width(rawX100));
fprintf('rawX98  : %d x %d\n',height(rawX98),width(rawX98));

fprintf('X100    : %d x %d\n',height(X100),width(X100));
fprintf('X98     : %d x %d\n',height(X98),width(X98));

fprintf('A100    : %d x %d\n',height(A100),width(A100));
fprintf('A98     : %d x %d\n',height(A98),width(A98));

fprintf('Y       : %d x %d\n',height(Y),width(Y));

fprintf('\nIPIP-100 domain sizes:\n');
fprintf('EXT = 20, EST = 20, AGR = 20, CSN = 20, OPN = 20\n');

fprintf('\nIPIP-98 domain sizes:\n');
fprintf('EXT = 20, EST = 20, AGR = 20, CSN = 20, OPN = 18\n');

fprintf('\nA100 variables:\n');
disp(string(A100.Properties.VariableNames)')

fprintf('\nA98 variables:\n');
disp(string(A98.Properties.VariableNames)')

fprintf('\nOutcome variables:\n');
disp(string(Y.Properties.VariableNames)')

%% ============================================================
% 16. Save CSV files
%% ============================================================

rawX100_file = fullfile(predictive_dir,'rawX100.csv');
rawX98_file = fullfile(predictive_dir,'rawX98.csv');

X100_file = fullfile(predictive_dir,'X100.csv');
X98_file = fullfile(predictive_dir,'X98.csv');

A100_file = fullfile(predictive_dir,'A100.csv');
A98_file = fullfile(predictive_dir,'A98.csv');

Y_file = fullfile(predictive_dir,'Y.csv');

writetable(rawX100,rawX100_file);
writetable(rawX98,rawX98_file);

writetable(X100,X100_file);
writetable(X98,X98_file);

writetable(A100,A100_file);
writetable(A98,A98_file);

writetable(Y,Y_file);

%% ============================================================
% Finished
%% ============================================================

fprintf('\n============================================================\n');
fprintf('FILES SAVED\n');
fprintf('============================================================\n');

fprintf('rawX100 : %s\n',rawX100_file);
fprintf('rawX98  : %s\n',rawX98_file);

fprintf('X100    : %s\n',X100_file);
fprintf('X98     : %s\n',X98_file);

fprintf('A100    : %s\n',A100_file);
fprintf('A98     : %s\n',A98_file);

fprintf('Y       : %s\n',Y_file);

fprintf('============================================================\n');