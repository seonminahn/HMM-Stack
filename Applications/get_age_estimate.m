function get_age_estimate(coreName, age_top, age_bottom)
%
% example: 
%     get_age_estimate('GeoB1032_LR04age', 1.42, 380)
%     get_age_estimate('GeoB1214_LR04age', 2.09, 450)
%

% Criteria of convergence
iterTol = 0.1;  %% Iteration stops if log likelihood does not increase larger than iterTol%
iterMax = 10;   %% Iteration stops if the number of iteration is larget than iterMax.


% Load a record and save the information of the record as a mat file
inputD = load(['data/' coreName '.txt']);
inputNAME{1} = ['../data/' coreName '.txt'];
beginNend = [1, size(inputD,1), age_top, age_bottom];
save(['data/' coreName '_info.mat'], 'inputNAME', 'beginNend')


% Get age estimates
cd code

contiIter = true;
bb = 0;
while contiIter == true
    clearvars -except forII age_bottom coreName contiIter iterTol iterMax

    run_for_communication(['age_estimate_' coreName], age_bottom, date, -1, 'step', 0);
    
    coreNAME_woL = ['age_estimate_' coreName]; 
    dataL = age_bottom; 
    inputN = 1; 
    r_date = date; 
    exactStart = -1; 
    target_est = 'step'; 
    dataNormalize = 0; 
    run_for_each_input;
    
    if bb ~= 0
        current_log_fE = log_fE;
        load(['../../Results/' date '/age_estimate_' coreName '_input1_iter' num2str(bb) '_log_fE.mat']);     
        pre_log_fE = log_fE;
%         disp(abs(current_log_fE-pre_log_fE)/pre_log_fE)
        if abs(current_log_fE-pre_log_fE)/pre_log_fE < iterTol/100
            contiIter = false;
        end
    end 
    
    if bb > iterMax
        contiIter = false;
    end
end


% Plot median estimates and age uncertainty. 
clearvars -except coreName
resultNS = ['age_estimate_' coreName '_input1_iter']; 
resultS = ['../../Results/' date '/age_estimate_' coreName '_input1_iter*'];
calc_median_ali_dp_fine_half; 
plot_align;
save(['../../Results/' date '/age_estimate_' coreName '_fullData'])

core_input = input_scaled_new(:,2);
core_median = median_ali;
core_upper95 = upper_95;
core_lower95 = lower_95;
core_ratio = ali_ratio;
save(['../../Results/' date '/age_estimate_' coreName '_summary'], 'core_input', 'core_median', 'core_upper95', 'core_lower95', 'core_ratio')

cd ..
