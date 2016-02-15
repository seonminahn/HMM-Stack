% function run_for_each_input(coreNAME_woL, dataL, inputN, r_date, exactStart, target_est, dataNormalize)
% clear all
% coreNAME_woL = 'LR04';
% dataL = 6000;
% inputN = 73;
% r_date = '05-May-2015';
% exactStart = 1;
% target_est = 'step';


%% script-based parallel test
% fid = fopen('script_parallel_test.txt', 'a+');
% fprintf(fid, 'Here is the inside of run_for_each_input(%d)\n', inputN);
% fclose(fid);

% target_est = 'step';
% target_est = 'linear';

% if exist('r_date', 'var') == 0
%     r_date = date;
% end
% if exist('exactStart', 'var') == 0
%     exactStart = 4;
% end
% if exist('target_est', 'var') == 0
%     target_est = 'linear';
% end
if exist('dataNormalize', 'var') == 0
    dataNormalize = true;
end

tic
disp(['input' num2str(inputN)])

coreNAME = coreNAME_woL;

resultsFolder = ['../../Results/', r_date, '/'];
updateFileN = [coreNAME '_iter'];
inputFileN  = [coreNAME '_input' num2str(inputN) '_iter'];

preD = dir([resultsFolder updateFileN '*']);
preDN = zeros(length(preD),1);
for i = 1 : length(preD)
    preDN(i) = sscanf(preD(i).name,  [updateFileN '%d']);
end

% disp(['Loading ' resultsFolder updateFileN num2str(max(preDN))])
load([resultsFolder updateFileN num2str(max(preDN)) '_updateD'])

    
if max(preDN) == 0        
    %% sequence_info
    ori1 = load(inputNAME{inputN});
    if strcmp(coreNAME_woL(1:3), 'sam')
%         total_target = load('sample_LR04.txt');
        ori2 = load(['../NewBenthicdata/' coreNAME_woL '/' coreNAME_woL '.txt']);
    else
        if strcmp(coreNAME_woL, 'LR04150')
            ori2 = load('LR04150');
        elseif strcmp(coreNAME_woL, 'LR04150f')
            ori2 = load('LR04150f');
        else
            ori2 = load('LR04_cib');
        end
        
    end
    ori2(:,2) = 0;
    
    begin1 = beginNend(inputN, 1);
    end1   = beginNend(inputN, 2);
    begin2 = beginNend(inputN, 3);
    end2   = beginNend(inputN, 4);
    
%     seq1 = ori1(min(find(ori1(:,2)>=begin1)):max(find(ori1(:,2)<=end1)),[1 3]);
%     seq2 = ori2(min(find(ori2(:,1)>=begin2)):max(find(ori2(:,1)<=end2)),:);
    seq1 = ori1(begin1:end1, [1 3]);
    seq2 = ori2(find(ori2(:,1)>=begin2,1,'first'):find(ori2(:,1)<=end2,1,'last'),:);
    
    if dataNormalize == true
        seq1(:,2) = seq1(:,2) - diff_from_grand_mean(inputN);
    end
    
    %% Preprocessingx
    pre_processing
    
    
else    
%     disp(['Loading ' resultsFolder inputFileN num2str(bb)])   
    load([resultsFolder inputFileN num2str(bb)])
    
    bb = bb + 1;

%     disp(['Loading ' resultsFolder updateFileN num2str(bb)])
    load([resultsFolder updateFileN num2str(bb) '_updateD'])

    if strcmp(target_est, 'step')
        for i = 1 : length(bin_mean_index)
            mu(i) = total_mu(bin_mean_index(i));
            sigma(i) = total_sigma(bin_mean_index(i));
        end
    elseif strcmp(target_est, 'linear')
        for i = 1 : length(bin_mean_index)
            mu(i) = bin_mean_index_p(i,2)*total_mu(bin_mean_index_p(i,1)) + bin_mean_index_n(i,2)*total_mu(bin_mean_index_n(i,1));
            sigma(i) = bin_mean_index_p(i,2)*bin_mean_index_p(i,2)*total_sigma(bin_mean_index_p(i,1)) + bin_mean_index_n(i,2)*bin_mean_index_n(i,2)*total_sigma(bin_mean_index_n(i,1));
        end
    else
        disp('SHOULD NOT REACH HERE')
    end
    mean_mu = mean(mu);
    mean_sigma = mean(sigma);
    
    mu = mu - meanShift;
    mean_mu = mean_mu - meanShift;
    
    log_grid = log(normpdf(delta_grid, ones(length(input_scaled_new),1)*mu', ones(length(input_scaled_new),1)*sigma'));
    
    log_input_gap_raw = zeros(length(input_scaled_new),1);
    for i=1:length(input_scaled_new)
        tempNormpdf = normpdf(input_scaled_new(i,2)-mean_tar, mean_mu, mean_sigma);
        log_input_gap_raw(i) = log(tempNormpdf);
    end
end

forward_point_m

% save('saveBeforeBack_sample')

% log_fE
back_sample_m

clear log_f

update_parameter



clear fhandle
runTime = toc;

mem = mem_use();
% disp(mem)

save([resultsFolder inputFileN num2str(bb+1)])
save([resultsFolder inputFileN num2str(bb+1) '_log_fE'], 'log_fE')
if strcmp(target_est, 'linear')
    save([resultsFolder inputFileN num2str(bb+1) '_updateD'], 'bin_mean_index', 'bin_mean_index_p', 'bin_mean_index_n', 'sample_seq', 'input_scaled_new')
else
    save([resultsFolder inputFileN num2str(bb+1) '_updateD'], 'bin_mean_index', 'sample_seq', 'input_scaled_new', 'meanShift')
end


% save([resultsFolder inputFileN num2str(bb+1) '_mem'], 'mem')

