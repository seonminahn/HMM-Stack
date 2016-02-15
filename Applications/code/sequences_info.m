function [inputNAME, beginNend] = sequences_info(coreNAME, dataL)

switch coreNAME
            
    case {'LR04_925'}
        inputNAME{1} = '../NewBenthicdata/925_LR04age.txt';
        beginNend(1,:) = [1 1687 0 5000];
        

    case {'LR04_MD052920'}
        inputNAME{1} = '../NewBenthicdata/NewRecords/MD05-2920.txt';
        beginNend(1,:) = [1 259 0 389];        
        
    case {'LR04', 'LR04_500', 'LR04_1000', 'LR04plusNew', 'LR04_2','LR04_5','LR04150', 'LR04150f'}
        
        oceanN = '../NewBenthicdata/LR04/';
        pList = dir([oceanN '*.txt']);
        
        beginD = 0;
        endD = dataL;
        count = 0;
        for i = 1 : length(pList)
            
            d = load([oceanN pList(i).name]);
%             seq_temp = d(min(find(d(:,2)>=begin1)):max(find(d(:,2)<=end1)),[1 3]);
            begin1 = find(d(:,2)>=beginD, 1, 'first');
            end1 = find(d(:,2)<=endD, 1, 'last');
            seq_temp = d(begin1:end1, [1 3]);
            
            clear seq
            [seq(:,1),label]=unique(seq_temp(:,1));
            seq(:,2)=seq_temp(label,2);
            
            if size(seq,1) ~= 0
                count = count + 1;
                inputNAME_temp{count,1} = [oceanN, pList(i).name];
                l(count,1) = length(seq);
                
                
%                 beginNend_temp(count, [1,3]) = floor(d(min(find(d(:,2)>=begin1)),2));
%                 beginNend_temp(count, [2,4]) = ceil(d(max(find(d(:,2)<=end1)),2)); 

                beginNend_temp(count, 1) = begin1;              % Column index of the first point
                beginNend_temp(count, 2) = end1;                % Column index of the last point
                beginNend_temp(count, 3) = floor(d(begin1,2));  % Age estimate of the first point
                beginNend_temp(count, 4) = ceil(d(end1,2));     % Age estimate of the last point

            end
        end
        
        [ll, li] = sort(l,1,'descend');
        for i = 1 : length(ll)
            inputNAME{i} = inputNAME_temp{li(i),1};
            beginNend(i,:) = beginNend_temp(li(i),:);
        end
        
%         if dataL == 6000
%             tmp1 = inputNAME;
%             tmp2 = beginNend;
%             
%             clear inputNAME beginNend
%             %         newInd = [1:16 37 17:27 29 28 30:36 39:70 72 73 71 38];
%             newInd = [1:27 29 28 30:69 71 72 70];
%             for i = 1 : length(newInd)
%                 inputNAME{i} = tmp1{newInd(i)};
%                 beginNend(i,:) = tmp2(newInd(i),:);
%             end
%         end
        
        if strcmp(coreNAME, 'LR04plusNew')
           [inputNAME_new, beginNend_new] = sequences_info_newRecords_whole(dataL);
            inputNAME = [inputNAME inputNAME_new];
            beginNend = [beginNend; beginNend_new];
        end
        
        if strcmp(coreNAME, 'LR04_2')
           inputNAME(1:end-2) = [];
           beginNend(1:end-2,:) = [];
        end
 
        if strcmp(coreNAME, 'LR04_5')
           inputNAME([1,7:end]) = [];
           beginNend([1,7:end],:) = [];
        end
        
    case {'sample_LR04', 'sample_offset_LR04'}
        
        oceanN = ['../NewBenthicdata/' coreNAME '/'];
        pList = dir([oceanN '*.txt']);
        
        beginD = 0;
        endD = dataL;
        count = 0;
        for i = 1 : length(pList)
            
            d = load([oceanN pList(i).name]);
%             seq_temp = d(min(find(d(:,2)>=begin1)):max(find(d(:,2)<=end1)),[1 3]);
            begin1 = find(d(:,2)>=beginD, 1, 'first');
            end1 = find(d(:,2)<=endD, 1, 'last');
            seq_temp = d(begin1:end1, [1 3]);
            
            clear seq
            [seq(:,1),label]=unique(seq_temp(:,1));
            seq(:,2)=seq_temp(label,2);
            
            if size(seq,1) ~= 0
                count = count + 1;
                inputNAME_temp{count,1} = [oceanN, pList(i).name];
                l(count,1) = length(seq);
                
                
%                 beginNend_temp(count, [1,3]) = floor(d(min(find(d(:,2)>=begin1)),2));
%                 beginNend_temp(count, [2,4]) = ceil(d(max(find(d(:,2)<=end1)),2)); 

                beginNend_temp(count, 1) = begin1;              % Column index of the first point
                beginNend_temp(count, 2) = end1;                % Column index of the last point
                beginNend_temp(count, 3) = floor(d(begin1,2));  % Age estimate of the first point
                beginNend_temp(count, 4) = ceil(d(end1,2));     % Age estimate of the last point

            end
        end
        
        [ll, li] = sort(l,1,'descend');
        for i = 1 : length(ll)
            inputNAME{i} = inputNAME_temp{li(i),1};
            beginNend(i,:) = beginNend_temp(li(i),:);
        end
        
        
        %% 32 Atlantic cores of LR04 (defined at sequenceLoad.m)
%         load('../NewBenthicdata/LR04/LR04_Atlantic.mat')
    case 'PA'
        inputNAME{1} = '../LR04cores_spec_corr/1012_LR04age.txt';
        inputNAME{2} = '../LR04cores_spec_corr/1020_LR04age.txt';
        inputNAME{3} = '../LR04cores_spec_corr/V21146_LR04age.txt';
        inputNAME{4} = '../LR04cores_spec_corr/1143_LR04age.txt';
        inputNAME{5} = '../LR04cores_spec_corr/V1928_LR04age.txt';
        inputNAME{6} = '../LR04cores_spec_corr/677_LR04age.txt';
        inputNAME{7} = '../LR04cores_spec_corr/RC13110_LR04age.txt';
        inputNAME{8} = '../LR04cores_spec_corr/846_LR04age.txt';
        inputNAME{9} = '../LR04cores_spec_corr/982_LR04age.txt';
        inputNAME{10} = '../LR04cores_spec_corr/1123_LR04age.txt';
        inputNAME{11} = '../LR04cores_spec_corr/1089_LR04age.txt';
        inputNAME{12} = '../LR04cores_spec_corr/980_LR04age.txt';
        inputNAME{13} = '../LR04cores_spec_corr/984_LR04age.txt';
        inputNAME{14} = '../LR04cores_spec_corr/658_LR04age.txt';
        inputNAME{15} = '../LR04cores_spec_corr/1090_LR04age.txt';
        inputNAME{16} = '../LR04cores_spec_corr/925_LR04age.txt';
        inputNAME{17} = '../LR04cores_spec_corr/927_LR04age.txt';
        inputNAME{18} = '../LR04cores_spec_corr/GeoB1041_LR04age.txt';
        inputNAME{19} = '../LR04cores_spec_corr/659_LR04age.txt';
        inputNAME{20} = '../LR04cores_spec_corr/664_LR04age.txt';
        
    case 'P'
        inputNAME{1} = '../LR04cores_spec_corr/1012_LR04age.txt';
        inputNAME{2} = '../LR04cores_spec_corr/1020_LR04age.txt';
        inputNAME{3} = '../LR04cores_spec_corr/V21146_LR04age.txt';
        inputNAME{4} = '../LR04cores_spec_corr/1143_LR04age.txt';
        inputNAME{5} = '../LR04cores_spec_corr/V1928_LR04age.txt';
        inputNAME{6} = '../LR04cores_spec_corr/677_LR04age.txt';
        inputNAME{7} = '../LR04cores_spec_corr/RC13110_LR04age.txt';
        inputNAME{8} = '../LR04cores_spec_corr/846_LR04age.txt';
        inputNAME{9} = '../LR04cores_spec_corr/982_LR04age.txt';
        inputNAME{10} = '../LR04cores_spec_corr/1123_LR04age.txt';
        
    case 'A'
        inputNAME{1} = '../LR04cores_spec_corr/1089_LR04age.txt';
        inputNAME{2} = '../LR04cores_spec_corr/980_LR04age.txt';
        inputNAME{3} = '../LR04cores_spec_corr/984_LR04age.txt';
        inputNAME{4} = '../LR04cores_spec_corr/658_LR04age.txt';
        inputNAME{5} = '../LR04cores_spec_corr/1090_LR04age.txt';
        inputNAME{6} = '../LR04cores_spec_corr/925_LR04age.txt';
        inputNAME{7} = '../LR04cores_spec_corr/927_LR04age.txt';
        inputNAME{8} = '../LR04cores_spec_corr/GeoB1041_LR04age.txt';
        inputNAME{9} = '../LR04cores_spec_corr/659_LR04age.txt';
        inputNAME{10} = '../LR04cores_spec_corr/664_LR04age.txt';
        
    case 'P2'
        inputNAME{1} = '../LR04cores_spec_corr/982_LR04age.txt';
        inputNAME{2} = '../LR04cores_spec_corr/1123_LR04age.txt';
        
    case 'P1'
        inputNAME{1} = '../LR04cores_spec_corr/982_LR04age.txt';

    case {'S1', 'S1_1'}
        for inputN = 1 : 20
           inputNAME{inputN} = ['../ProfileHMM_generative_model/samples/input_woID_0to150_', num2str(inputN), '.txt']; 
        end
    case 'S2'
        for inputN = 1 : 40
           inputNAME{inputN} = ['../ProfileHMM_generative_model/samples/input_woID_0to150_', num2str(inputN), '.txt']; 
        end
    case {'S3', 'S3_1'}
        for inputN = 1 : 20
           inputNAME{inputN} = ['../ProfileHMM_generative_model/samples/input_woID_0to150_', num2str(inputN+40), '.txt']; 
        end    
        
    case 'S300_0'
        for inputN = 1 : 40
           inputNAME{inputN} = ['../ProfileHMM_generative_model/samples/input_woID_0to300_int0_', num2str(inputN), '.txt'];
        end
    case 'S300_1'
        for inputN = 1 : 40
           inputNAME{inputN} = ['../ProfileHMM_generative_model/samples/input_woID_0to300_int1_', num2str(inputN), '.txt'];
        end
    case 'S300_2'
        for inputN = 1 : 40
           inputNAME{inputN} = ['../ProfileHMM_generative_model/samples/input_woID_0to300_int2_', num2str(inputN), '.txt'];
        end
    case 'S300_4'
        for inputN = 1 : 40
           inputNAME{inputN} = ['../ProfileHMM_generative_model/samples/input_woID_0to300_int4_', num2str(inputN), '.txt'];
        end   
        
    case {'S300_5', 'S300_15', 'S300_25', 'S300_35', 'S300_45', 'S300_55'}
        for inputN = 1 : 20
           inputNAME{inputN} = ['../ProfileHMM_generative_model/samples/input_woID_0to300_int0_', num2str(inputN), '.txt'];
        end
    case {'S300_6', 'S300_16', 'S300_26', 'S300_36', 'S300_46', 'S300_56'}
        for inputN = 1 : 20
           inputNAME{inputN} = ['../ProfileHMM_generative_model/samples/input_woID_0to300_int1_', num2str(inputN), '.txt'];
        end
    case {'S300_7', 'S300_17', 'S300_27', 'S300_37', 'S300_47', 'S300_57'}
        for inputN = 1 : 20
           inputNAME{inputN} = ['../ProfileHMM_generative_model/samples/input_woID_0to300_int2_', num2str(inputN), '.txt'];
        end
    case {'S300_8', 'S300_18', 'S300_28', 'S300_38', 'S300_48', 'S300_58'}
        for inputN = 1 : 20
           inputNAME{inputN} = ['../ProfileHMM_generative_model/samples/input_woID_0to300_int4_', num2str(inputN), '.txt'];
        end   
end

if exist('beginNend', 'var') == 0
    num_seq = size(inputNAME,2);
    
    %     beginNendseed = [0 dataL 0 dataL];
    beginD = 0;
    endD = dataL;
    beginNend = zeros(num_seq, 4);
    for i = 1 : num_seq
        d = load(inputNAME{i});
        begin1 = find(d(:,2)>=beginD, 1, 'first');
        end1 = find(d(:,2)<=endD, 1, 'last');
        
%         beginNend(i,:) = beginNendseed;
        beginNend(i, 1) = begin1;              % Column index of the first point
        beginNend(i, 2) = end1;                % Column index of the last point
        beginNend(i, 3) = floor(d(begin1,2));  % Age estimate of the first point
        beginNend(i, 4) = ceil(d(end1,2));
    end
end

% inputNAME{5} = 'LR04cores_spec_corr/1012_LR04age.txt';
% inputNAME{6} = 'LR04cores_spec_corr/RC13110_LR04age.txt';
% inputNAME{7} = 'LR04cores_spec_corr/V21146_LR04age.txt';
% inputNAME{8} = 'LR04cores_spec_corr/V1928_LR04age.txt';
% inputNAME{1} = 'LR04cores_spec_corr/1020_LR04age.txt'
% inputNAME{2} = 'LR04cores_spec_corr/1146_LR04age.txt'