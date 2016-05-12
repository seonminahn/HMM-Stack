function sum_log_fE = run_for_communication(coreNAME_woL, dataL, r_date, exactStart, target_est, dataNormalize)

tic
disp(' ')

coreNAME = coreNAME_woL;

resultsFolder = ['../../Results/', r_date, '/'];
if exist(resultsFolder, 'dir') ~= 7
    mkdir(resultsFolder)
end
updateFileN = [coreNAME '_iter'];
preD = dir([resultsFolder updateFileN '*']);


%% Load data file or Initialize parameters
if isempty(preD)    
    
    if exist('discreteRatio', 'var') == 0
        discreteRatio = true;
%         discreteRatio = false;
    end
    
    if exist('logComp', 'var') == 0
        logComp = true;
%         logComp = false;
    end
    
%     dataL = 500;
    if exactStart == -1
        inputInfo = ['../data/' strrep(coreNAME, 'age_estimate_', '') '_info.mat'];
        load(inputInfo)
    else
        if strcmp(coreNAME_woL, 'newStack') == 1
            inputInfo = '../recordSummary_info.mat';
            load(inputInfo)
        else
            [inputNAME, beginNend] = sequences_info(coreNAME_woL, dataL);
        end
    end
    num_seq = size(inputNAME,2);
        
    beginNendseed = [0 dataL 0 dataL];
%     beginNend = zeros(num_seq, 4);
%     for i = 1 : num_seq
%         beginNend(i,:) = beginNendseed;
%     end
    
    if strcmp(coreNAME_woL(1:3), 'sam')
%         total_target = load('sample_LR04.txt');
        total_target = load(['../NewBenthicdata/' coreNAME_woL '/' coreNAME_woL '.txt']);
    else
        if strcmp(coreNAME_woL, 'LR04150')
            total_target = load('LR04150');
        elseif strcmp(coreNAME_woL, 'LR04150f')
            total_target = load('LR04150f');
        else
            total_target = load('LR04_cib');
        end
    end
%     total_target = total_target( min(find(total_target(:,1)>=beginNendseed(3))) :max(find(total_target(:,1)<=beginNendseed(4))),:);
    total_target = total_target( find(total_target(:,1)>=beginNendseed(3),1,'first') : find(total_target(:,1)<=beginNendseed(4),1,'last'),:);
    ori_mu = total_target(:,2);
    total_target(:,2) = 0;
    
    
    
    %% Compute the mean of entire inputs
    % Scale inputs to have the same average -- Didn't consider the range of
    % input
%     seq1avg = zeros(num_seq,1);
%     for i = 1 : num_seq
%         ori1 = load(inputNAME{i});
% %         seq1 = ori1(min(find(ori1(:,2)>=beginNend(i,1))):max(find(ori1(:,2)<=beginNend(i,2))),[1 3]);
%         seq1 = ori1(beginNend(i,1):beginNend(i,2),[1 3]);
%         seq1avg(i) = sum(seq1(:,2))/length(seq1(:,2));
%     end
%     total_seq1avg = mean(seq1avg);
%     diff_from_grand_mean = seq1avg - total_seq1avg;

    % Scale inputs to have the same average with LR04 
    seq1avg = zeros(num_seq,1);
    seq2avg = zeros(num_seq,1);
    for i = 1 : num_seq
        begin1 = beginNend(i, 1);
        end1   = beginNend(i, 2);
        begin2 = beginNend(i, 3);
        end2   = beginNend(i, 4);
        
        ori1 = load(inputNAME{i});
        seq1 = ori1(begin1:end1, [1 3]);
        seq2 = ori_mu(find(total_target(:,1)>=begin2,1,'first'):find(total_target(:,1)<=end2,1,'last'));
        
        seq1avg(i) = sum(seq1(:,2))/length(seq1(:,2));
        seq2avg(i) = sum(seq2)/length(seq2);
    end
    total_seq1avg = mean(ori_mu);
        
    
    if dataNormalize == true
        diff_from_grand_mean = seq1avg - seq2avg;
        ori_mu_normalized = ori_mu - (mean(ori_mu) - total_seq1avg);
    else
        diff_from_grand_mean = 0*seq2avg;
        ori_mu_normalized = ori_mu;
    end
    
    
    %% Set IC for the new stack
    if exactStart == 1 % IC = LR04
        total_mu = ori_mu_normalized;             
        
    elseif exactStart == 10
        total_mu = load('sample_0_500_linear_model1.txt');
        
    elseif exactStart == 20
        total_mu = load('sample_0_500_linear_model2.txt');
        
    elseif exactStart == 980
        total_mu = load('IC_980.txt');
    elseif exactStart == 1148
        total_mu = load('IC_1148.txt');
        
    elseif exactStart == 5001
        total_mu = load('IC_sample1.txt');
    elseif exactStart == 50040
        total_mu = load('IC_sample40.txt');      
        
    elseif exactStart == 11
        total_mu = load('LR04_cib_linear_model_0_500.txt');
        
%     elseif exactStart > 100
%         total_mu = cos(2*pi/(exactStart/10)*total_target(:,1))+mean(ori_mu_normalized);
% 
%     elseif exactStart > 10
%         total_mu = sin(2*pi/exactStart*total_target(:,1))+mean(ori_mu_normalized);
%         
%         
%     elseif exactStart > 1 % IC = input having the highest resolution
%         if exactStart == 2
%             inputIC = load('../LR04cores_spec_corr/1089_LR04age.txt');
%         elseif exactStart == 3
%             inputIC = load('../LR04cores_spec_corr/980_LR04age.txt');
%         elseif exactStart == 4
%             inputIC = load(inputNAME{1});
%         elseif exactStart == 5
%             inputIC = load(inputNAME{2});
%         end
%         
% %         inputIC_mu = inputIC(min(find(inputIC(:,2)>=beginNendseed(3))):max(find(inputIC(:,2)<=beginNendseed(4))),[2,3]);
%         inputIC_mu = inputIC(find(inputIC(:,2)>=beginNendseed(3),1,'first'):find(inputIC(:,2)<=beginNendseed(4),1,'last'),[2,3]);
%         if dataNormalize == true
%             inputIC_mu(:,2) = inputIC_mu(:,2) - (mean(inputIC_mu(:,2))-total_seq1avg);
%         end
%         
%         total_mu = zeros(size(total_target,1), 1);
%         for targetI = 1 : size(total_target,1)
%             [~,tempI] = min(abs(inputIC_mu(:,1) - total_target(targetI,1)));
%             total_mu(targetI) = inputIC_mu(tempI,2);
%         end
%         %     figure
%         %     hold on
%         %     plot(total_target(:,1), ori_mu, 'k*-')
%         %     plot(total_target(:,1),total_mu, 'bo-')
%         %     plot(input1089_mu(:,1), input1089_mu(:,2), 'r:*')  
        
    else % IC = mean(LR04)
        total_mu = 0*ori_mu_normalized + mean(ori_mu_normalized);
    end
    
    total_sigma = 0*ori_mu_normalized + 0.25;
    
    if exactStart == -1
        newstack = load('../../Prob_stack.txt');
        total_mu = newstack(:,2);
        total_sigma = newstack(:,3);
    end
       
    totalMU(:,1) = total_mu;
    totalSIGMA(:,1) = total_sigma;
    
    total_meanShift = zeros(num_seq,1);
    
    bb = 0;
    sum_log_fE_old = 0;
    totalT = 1;
    save([resultsFolder updateFileN num2str(bb)])
    save([resultsFolder updateFileN num2str(bb) '_updateD'], 'total_target', 'total_mu', 'total_sigma', 'bb', 'diff_from_grand_mean', 'beginNend', 'inputNAME', 'dataNormalize', 'discreteRatio', 'logComp', 'total_meanShift')
    disp(['Initialization is saved at ' updateFileN '0.mat'])
    
    sum_log_fE = 0;
else 
    
    preDN = zeros(length(preD),1);
    for i = 1 : length(preD)
        preDN(i) = sscanf(preD(i).name,  [updateFileN '%d']);
    end

%     disp(['Loading ' resultsFolder updateFileN num2str(max(preDN))])
    load([resultsFolder updateFileN num2str(max(preDN))])

    bb = bb + 1;
    
    count_completed_core = 0;
    for i = 1 : num_seq
        fN = [resultsFolder coreNAME '_input' num2str(i) '_iter' num2str(bb) '_updateD.mat'];
        if exist(fN, 'file')
            count_completed_core = count_completed_core + 1;
            completed_cores(count_completed_core) = i;
        else
            disp(['No result from ' num2str(i)])
        end
    end

    total_log_fE = zeros(num_seq,1);
    for i = completed_cores
        load([resultsFolder coreNAME '_input' num2str(i) '_iter' num2str(bb) '_log_fE'])
        total_log_fE(i) = log_fE;
    end
    
    sum_log_fE = sum(total_log_fE);
    diff_log_fE = sum_log_fE - sum_log_fE_old;
    sum_log_fE_old = sum_log_fE;
    record_log_fE(bb+1,1) = sum_log_fE;
    fprintf('Total log likelihood in iteration %d is %f\n',bb, sum_log_fE);
    fprintf('The increased amount of the total log likelihood is %f\n', diff_log_fE);
    
    if exactStart ~= -1
        disp('Updating mu and sigma')
        if diff_log_fE < sum_log_fE*0.0001
            totalT = 0;
            fprintf('Likelihood Converged\n')
            update_parameter_musigma
        else
            totalT = 1;
            
            update_parameter_musigma
        end
    end
    
    
    totalMU(:,bb+1) = total_mu;
    totalSIGMA(:,bb+1) = total_sigma;
    %totalAligned_record{bb+1} = aligned_record;
        
    save([resultsFolder updateFileN num2str(bb)])
    save([resultsFolder updateFileN num2str(bb) '_updateD'], 'total_mu', 'total_sigma', 'bb', 'dataNormalize', 'discreteRatio')
end


disp(['Communication is done for iter = ' num2str(bb)])
toc
disp(' ')

if strcmp(coreNAME_woL, 'newStack') == 1
total_record(:,1) = total_target(:,1);
total_record(:,2) = total_mu;
total_record(:,3) = total_sigma;
total_record(:,4) = total_mu+1.96*total_sigma;
total_record(:,5) = total_mu-1.96*total_sigma;

save([resultsFolder updateFileN num2str(bb) '.txt'], 'total_record', '-ascii')
end

