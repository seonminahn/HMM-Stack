# HMM-stack

HMM stack (HMM_stack.txt): Probabilistic stack constructed from 180 benthic *δ<sup>18</sup>O* records

HMM_LR04 stack (HMM_LR04_stack.txt): Probabilistic stack constructed from the LR04 cores

The stack files contain five columns: Age, *δ<sup>18</sup>O* value [‰], standard deviation of *δ<sup>18</sup>O* value, upper bound of 95% interval, and lower bound of 95% interval. 

## Applications

The Application folder includes three applications: construction of a probabilistic stack from benthic *δ<sup>18</sup>O* records, Age estimation of a benthic *δ<sup>18</sup>O* record, and lead/lag analysis between two events observed from different cores. 

All codes are written in MATLAB and located under the code folder. 

All *δ<sup>18</sup>O* record files should be located under the data folder. Each record should include three columns: depth, age, and data value. When age estimates are unavailable, you can leave them as NaN.

### Probabilistic stack

The MATLAB code 'construct_hmm_stack' constructs a probabilistic stack using cores listed in the summary file. To construct a probabilistic stack with your own benthic *δ<sup>18</sup>O* records, you should modify the 'recordSummary.txt' file, located under the data folder. Each line of this summary file includes three columns: the file name of the core, age estimates of the top and the bottom of the core. Then, you can run the MATLAB code 'construct_hmm_stack' with the summary file as follows. 

    construct_hmm_stack('recordSummary.txt')

This code generates 'Yourstack.txt', which is a new probabilistic stack constructed from the cores listed in 'recordSummary.txt'. It follows the same format as HMMstack.txt. In addition to the stack file, the code generates following files while running the codes. 

* YourStack_iterN.mat, YourStack_iterN_updateD.mat: data files saved after the N<sup>th</sup> iteration
* YourStack_inputM_iterN.mat, YourStack_inputM_iterN_updateD.mat: data files save after aligning the M<sup>th</sup> benthic *δ<sup>18</sup>O* record to the stack generated after the N<sup>th</sup> iteration

Do not delete any files while running codes. These files are required to update each iteration. 

### Age estimation

The MATLAB code 'get_age_estimate' finds age estimates of a benthic *δ<sup>18</sup>O* record by aligning the record with the HMM stack. To run the code, you should provide the file name of the core and age estimates of the top and the bottom of the core as follow.

    get_age_estimate('coreName', age_top, age_bottom)

The code generates 'coreName_HMMstack.mat' which includes the following variables:

* core_input: benthic *δ<sup>18</sup>O* values of the core
* core_median: median age estimates of the core 
* core_upper95: upper bounds 95% interval for age estimates
* core_lower95: lower bounds of 95% interval for age estimates
* core_ratio: relative accumulation ratio to the HMM stack based on the median age estimates


### Lead/Lag analysis

The MATLAB code 'analyze_lead_lag' provides the probability of one point in a record coreNAME1 leading to another point in a record coreNAME2. Before running this code, age estimates of coreNAME1 and coreNAME2 should be estimated using 'get_age_estimate'. To run this code, you should specify two points you want to analyze using their depths as follow. 

    analyze_lead_lag('coreName1', depth1, 'coreName2', depth2)


## License

The HMM-stack algorithm is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation.
