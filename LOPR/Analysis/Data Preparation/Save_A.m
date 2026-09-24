%% Save LOOPR latent A matrix

clear

clc

root_dir = "/Users/avinandanroy/Desktop/Research/DDE/DDE_Fitting_Pshycology";

load(fullfile(root_dir,"LOPR","Analysis","PostProcessed_Data","LOOPR_PostProcessed.mat"),"LOOPR");

A1 = LOOPR.A1;

A2 = LOOPR.A2;

A = [A1 A2];

A1_names = "A1_" + string(1:size(A1,2));

A2_names = "A2_" + string(1:size(A2,2));

A = array2table(A,'VariableNames',[A1_names A2_names]);

out_dir = fullfile(root_dir,"LOPR","Analysis","Predictive_Data");

if ~exist(out_dir,'dir')

    mkdir(out_dir);

end

writetable(A,fullfile(out_dir,"A.csv"));

fprintf('Saved A.csv\n');

fprintf('Dimensions: %d x %d\n',height(A),width(A));

fprintf('A1 factors: %d\n',size(A1,2));

fprintf('A2 factors: %d\n',size(A2,2));