function [anchor_table, anchor_question_table, top_tmp_values] = ...
    find_key_item_IPIPFMM(B1_CSP, top_N)

% FIND_KEY_ITEM_IPIPFMM
%
% Finds anchor (key) questionnaire items for each active first-layer factor
% in the IPIP-FFM personality data.
%
% INPUTS
% -------
% B1_CSP : J x (K1+1) loading matrix (first column = intercept)
% top_N  : Number of anchor items to return (default = 3)
%
% OUTPUTS
% --------
% anchor_table           : Table of anchor variable names
% anchor_question_table  : Table of corresponding questionnaire text
% top_tmp_values         : Anchor strength values

if nargin < 2
    top_N = 3;
end

%% -------------------------------------------------------------
%% Variable names
%% -------------------------------------------------------------

survey_vars = [
    compose("EXT%d",1:10), ...
    compose("EST%d",1:10), ...
    compose("AGR%d",1:10), ...
    compose("CSN%d",1:10), ...
    compose("OPN%d",1:10)];

var_names = cellstr(survey_vars);

question_text = [
    "I am the life of the party."
    "I don't talk a lot."
    "I feel comfortable around people."
    "I keep in the background."
    "I start conversations."
    "I have little to say."
    "I talk to a lot of different people at parties."
    "I don't like to draw attention to myself."
    "I don't mind being the center of attention."
    "I am quiet around strangers."

    "I get stressed out easily."
    "I am relaxed most of the time."
    "I worry about things."
    "I seldom feel blue."
    "I am easily disturbed."
    "I get upset easily."
    "I change my mood a lot."
    "I have frequent mood swings."
    "I get irritated easily."
    "I often feel blue."

    "I feel little concern for others."
    "I am interested in people."
    "I insult people."
    "I sympathize with others' feelings."
    "I am not interested in other people's problems."
    "I have a soft heart."
    "I am not really interested in others."
    "I take time out for others."
    "I feel others' emotions."
    "I make people feel at ease."

    "I am always prepared."
    "I leave my belongings around."
    "I pay attention to details."
    "I make a mess of things."
    "I get chores done right away."
    "I often forget to put things back in their proper place."
    "I like order."
    "I shirk my duties."
    "I follow a schedule."
    "I am exacting in my work."

    "I have a rich vocabulary."
    "I have difficulty understanding abstract ideas."
    "I have a vivid imagination."
    "I am not interested in abstract ideas."
    "I have excellent ideas."
    "I do not have a good imagination."
    "I am quick to understand things."
    "I use difficult words."
    "I spend time reflecting on things."
    "I am full of ideas."
    ];

var_mapping = table( ...
    var_names', ...
    question_text, ...
    'VariableNames', ...
    {'Variable','Question'});

%% -------------------------------------------------------------
%% Active factors
%% -------------------------------------------------------------

J = size(B1_CSP,1);

active_factor_index = ...
    find(any(B1_CSP(:,2:end)~=0,1));

active_cols_B1 = active_factor_index + 1;

K1_active = length(active_cols_B1);

top_items_B1 = strings(K1_active,top_N);
top_tmp_values = zeros(K1_active,top_N);

%% -------------------------------------------------------------
%% Find anchor items
%% -------------------------------------------------------------

for kk = 1:K1_active

    k_col = active_cols_B1(kk);

    other_cols = setdiff(active_cols_B1,k_col);

    tmp = zeros(J,1);

    for j = 1:J

        if isempty(other_cols)

            tmp(j) = max(0,B1_CSP(j,k_col));

        else

            diff_load = ...
                B1_CSP(j,k_col) - ...
                B1_CSP(j,other_cols);

            tmp(j) = max(0,min(diff_load));

        end

    end

    [sorted_tmp,sort_idx] = sort(tmp,'descend');

    keep = find(sorted_tmp>0);

    n_keep = min(top_N,length(keep));

    if n_keep>0

        idx = sort_idx(1:n_keep);

        top_items_B1(kk,1:n_keep) = ...
            string(var_names(idx));

        top_tmp_values(kk,1:n_keep) = ...
            sorted_tmp(1:n_keep);

    end

end

%% -------------------------------------------------------------
%% Build tables
%% -------------------------------------------------------------

keep_factor = any(top_tmp_values>0,2);

factor_names = ...
    "A1_" + string(active_factor_index(keep_factor));

anchor_table = array2table( ...
    top_items_B1(keep_factor,:)', ...
    'VariableNames',cellstr(factor_names), ...
    'RowNames',cellstr("Key Item " + string(1:top_N)));

anchor_question_table = anchor_table;

for c = 1:width(anchor_question_table)

    for r = 1:height(anchor_question_table)

        item = anchor_question_table{r,c};

        if strlength(item)>0

            idx = strcmp(var_mapping.Variable,item);

            anchor_question_table{r,c} = ...
                var_mapping.Question(idx);

        end

    end

end

end