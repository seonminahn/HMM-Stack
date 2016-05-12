# Probabilistic stacks

Prob-stack (Prob_stack.txt): Probabilistic stack constructed from 180 benthic *δ<sup>18</sup>O* records. (Detailed information about the records is given in the metatdata table of [1].)

Prob-LR04-stack (Prob_LR04_stack.txt): Probabilistic stack constructed from the LR04 cores.

Both stack files contain five columns: Age, *δ<sup>18</sup>O* value [‰], standard deviation of *δ<sup>18</sup>O* value, upper bound of 95% interval, and lower bound of 95% interval. 

## Applications

The Application folder includes three applications: construction of a probabilistic stack from benthic *δ<sup>18</sup>O* records, age estimation of a benthic *δ<sup>18</sup>O* record, and lead/lag analysis between two events observed from different cores. 

All codes are written in the MATLAB language and located in the code folder. 

All *δ<sup>18</sup>O* record files should be located in the data folder. To run this program, each record should include three columns: depth, age, and data value. When age estimates are unavailable, you can leave them as NaN.

To run these applications, please download the Application folder and the Prob-stack; for the program to run, they must remain in the same folder. 

### Probabilistic stack

The MATLAB code 'construct_hmm_stack' constructs a probabilistic stack using cores listed in the summary file. To construct a probabilistic stack with your own benthic *δ<sup>18</sup>O* records, you should modify the 'recordSummary.txt' file, located in the data folder. Each line of this summary file includes three columns: the file name of the core, the age estimate of the top of the core, and the age estimate of the bottom of the core. It is enough to input rough estimates of the ages. After modifying the file, you can run the MATLAB code 'construct_hmm_stack' with the summary file as follows: 

    construct_hmm_stack('recordSummary.txt')

This code generates 'newStack_iterN.txt', which is a new probabilistic stack constructed from the cores listed in 'recordSummary.txt'. It follows the same format as Prob_stack.txt. In addition to the stack file, while running, the code generates the following files: 

* newStack_iterN.mat, newStack_iterN_updateD.mat (data files that are saved after the N<sup>th</sup> iteration)
* newStack_inputM_iterN.mat, newStack_inputM_iterN_updateD.mat (data files that are saved after aligning the M<sup>th</sup> benthic *δ<sup>18</sup>O* record to the stack generated after the N<sup>th</sup> iteration)

Do not delete any files while running these codes. These files are required to update each iteration. 

To construct a valid probabilistic stack, at least two different *δ<sup>18</sup>O* valued samples should be assigned to each point in the stack. Otherwise, the variance term becomes too small to construct a valid emission model, especially with a small number of records. To avoid this, the algorithm manually assigns a reasonably large value when a variance term is too small and prints an error message. Such cases did not occur when we constructed the Prob-stack and the Prob-LR04-stack because we used sufficient records. 

Note that a construction of the stack requires extensive use of high performance computing resources. Constructing the Prob-stack took about one month with 180 computing nodes as it takes over 100 hours for a high-resolution record to complete one iteration. It took about 30 minutes to construct a stack on a desktop (2.7 GHz Intel Core i% processor & 16GB memory) using the three example records as provided in 'recordSummary.txt'.

### Age estimation

The MATLAB code 'get_age_estimate' finds age estimates of a benthic *δ<sup>18</sup>O* record by aligning the record with the Prob-stack. To run the code, you should provide the file name of the core and age estimates of the top and the bottom of the core as follows:

    get_age_estimate(coreName, age_top, age_bottom)

The code generates 'age_estimate_coreName_summary.mat' which includes the following variables:

* core_input (benthic *δ<sup>18</sup>O* values of the core)
* core_median (median age estimates of the core) 
* core_upper95 (upper bounds of 95% interval for age estimates)
* core_lower95 (lower bounds of 95% interval for age estimates)
* core_ratio (relative accumulation ratio to the Prob-stack based on the median age estimates)


### Lead/Lag analysis

The MATLAB code 'analyze_lead_lag' provides the probability of one point in a record "coreNAME1" leading to another point in a record "coreNAME2". Before running this code, the age estimates of "coreNAME1" and "coreNAME2" should be estimated using 'get_age_estimate'. To run this code, you should specify two points you want to analyze using their depths as follows:

    analyze_lead_lag(coreName1, depth1, coreName2, depth2)


## License

The HMM-stack algorithm is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation.

The following citation should be included in any research reports, papers, or publications based on the Prob-stack and the HMM-stack algorithm. 

[1] A probabilistic Pliocene-Pleistocene stack of benthic *δ<sup>18</sup>O* using a profile hidden Markov model
