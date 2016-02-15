clearvars
dataL = 6000
oceanN = '../NewBenthicdata/LR04/';
pList = dir([oceanN '*.txt']);

begin1 = 0;
end1 = dataL;
count = 0;
for ii = 1 : length(pList)
    d = load([oceanN pList(ii).name]);
    seq1 = d(min(find(d(:,2)>=begin1)):max(find(d(:,2)<=end1)),[1 3]);
    
    if size(seq1,1) ~= 0
        clear input
        
        count = count + 1;
        inputNAME_temp{count,1} = [oceanN, pList(ii).name];
        
        [input(:,1),label]=unique(seq1(:,1));
        input(:,2)=seq1(label,2);
        
        l(count,1) = length(input);
        im(count,1) = round(d(min(find(d(:,2)>=begin1)),2));
        iM(count,1) = round(d(max(find(d(:,2)<=end1)),2));
        marginBE(count,1) = 0; %round((iM(count)-im(count))*0.05);
        beginNend_temp(count,[1 3]) = max(0, im(count)-marginBE(count));
        beginNend_temp(count,[2 4]) = iM(count)+marginBE(count);
        
        tB = max(0, im(count)-marginBE(count));
        tE = iM(count)+marginBE(count);
        input_scaled_temp{count} = tB+(input(:,1)-input(1,1))*(tE-tB)/(input(end,1)-input(1,1));
        sortedDiff = sort(diff(input_scaled_temp{count}));
        targetDistance_temp(count) = max(sortedDiff(ceil(length(sortedDiff)*0.10))/4, 1/10);
        discrete_target_temp{count} = tB : targetDistance_temp(count) : tE;
        discrete_target_l_temp(count) = length(discrete_target_temp{count});
        if discrete_target_l_temp(count) == 0
            disp('here')
        end
        
        discrete_target = discrete_target_temp{count};
        input_scaled = input_scaled_temp{count};
        input_bin_index = zeros(length(input),1);
        ss =1;
        for i = 1: length(discrete_target)
            while input_scaled(ss,1)<=discrete_target(i)
                input_bin_index(ss)=i;
                ss = ss+1;
                if ss>length(input)
                    break
                end
            end
        end
        input_bin_index(ss:end)=length(discrete_target);
        
        
        %record the number of bin distance betwen input point with previous point
        %set the distance for the first element as 0
        dis_bin_input = zeros(length(input),1);
        for i=2:length(dis_bin_input)
            dis_bin_input(i) = input_bin_index(i) - input_bin_index(i-1);
        end
        maxbin_temp_whole = sort(dis_bin_input, 'descend');
        maxbin_temp(count,1) = maxbin_temp_whole(ceil(length(maxbin_temp_whole)*0.005));
        maxbin_temp_2(count,1) = max(maxbin_temp_whole);
        
        log_f = cell(length(input_scaled),length(discrete_target));
        for ilog = 1 : size(log_f,1)
            for jlog = 1 : size(log_f,2)
                log_f{ilog,jlog} = zeros(4*dis_bin_input(ilog),1);
            end
        end
        whos_log_f = whos('log_f');
        size_log_f_temp(count,1) = whos_log_f.bytes;
        
    end
end

[ll, li] = sort(l,1,'descend');
for i = 1 : length(ll)
    inputNAME{i} = inputNAME_temp{li(i),1};
    beginNend(i,:) = beginNend_temp(li(i),:);
    discrete_target_l(i,1) = discrete_target_l_temp(li(i));
    targetDistance(i,1) = targetDistance_temp(li(i));
    input_scaled_new{i,1} = input_scaled_temp{li(i)};
    discrete_target_new{i,1} = discrete_target_temp{li(i)};
    im_new(i,1) = im(li(i));
    iM_new(i,1) = iM(li(i));
    numPerOne(i,1) = discrete_target_l(i)/(iM_new(i) - im_new(i));
    maxbin(i,1) = maxbin_temp(li(i));
    maxbin_2(i,1) = maxbin_temp_2(li(i));
    size_log_f(i,1) = size_log_f_temp(li(i))/1024/1024/1024;
end


for i = 1 : length(ll)
    a(i,1) = max(maxbin_2(i), maxbin(i)*4)*discrete_target_l(i)*length(input_scaled_new{i})*8/1024/1024/1024;
    b(i,1) = max(maxbin(i)*4)*discrete_target_l(i)*length(input_scaled_new{i})*8/1024/1024/1024;
    o(i,1) = maxbin_2(i)*4*discrete_target_l(i)*length(input_scaled_new{i})*8/1024/1024/1024;
end

c = a - b;
[cc ci] = sort(c, 'descend');

% %%
% inputN = 1;
% thisInput = input_scaled{inputN,1};
% diffI = diff(thisInput);
% sdiffI = sort(diffI);
%

if dataL >= 6000
    tmp1 = inputNAME;
    tmp2 = beginNend;
    
    clear inputNAME beginNend
    % newInd = [1:16 24 17:23 25:27 29 28 30:37 39:70 72 73 71 38];
    newInd = [1:27 29 28 30:69 71 72 70];
    
    for i = 1 : length(newInd)
        inputNAME{i} = tmp1{newInd(i)};
        beginNend(i,:) = tmp2(newInd(i),:);
    end
end




