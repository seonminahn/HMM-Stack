%%
clear
coreNAME_woL = 'LR04_925'; dataL = 5000; inputN = 1; r_date = '30-Nov-2015'; exactStart = -1; target_est = 'step'; dataNormalize = 0; run_for_each_input;
run_for_communication('LR04_925', 5000, '30-Nov-2015', -1, 'step', 0);

clear 
coreNAME_woL = 'LR04_925'; dataL = 5000; inputN = 1; r_date = '30-Nov-2015'; exactStart = -1; target_est = 'step'; dataNormalize = 0; run_for_each_input;
run_for_communication('LR04_925', 5000, '30-Nov-2015', -1, 'step', 0);


%%
for forII = 1 : 10
    forII
clearvars except forII 
run_for_communication('LR04_MD052920', 500, 'Dec-8-2015', -1, 'step', 0);
coreNAME_woL = 'LR04_MD052920'; dataL = 500; inputN = 1; r_date = 'Dec-8-2015'; exactStart = -1; target_est = 'step'; dataNormalize = 0; run_for_each_input;
end

%%
resultNS = 'LR04_MD052920_500_IC-1_step_norm0_input1_bb'; 
resultS = '../../ProfileHMM_results/Dec-8-2015/LR04_MD052920_500_IC-1_step_norm0_input1_bb*';
calc_median_ali_dp_fine_half