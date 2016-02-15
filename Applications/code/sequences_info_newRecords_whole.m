function [inputNAME, beginNend] = sequences_info_newRecords_whole(dataL)
% clear all

% dataL = 1000;
testM = false;

% New records list
oceanN = '../NewBenthicdata/NewRecords/';
f_list_id = '../NewBenthicdata/NewRecords_Name.txt';
fid = fopen(f_list_id);

% Age estimate of new records
age_estimate = load('../NewBenthicdata/NewRecords_Age.txt');

% figure;

count = 0;
inputNAME_temp = cell(length(age_estimate),1);
beginNend_temp = zeros(length(age_estimate),4);
l = zeros(length(age_estimate),1);

while ~feof(fid)
    tline = fgetl(fid)
    count = count + 1;
    inputNAME_temp{count,1} = [oceanN tline '.txt'];
    
    d = load(inputNAME_temp{count,1});
    beginNend_temp(count, 1) = 1;
    beginNend_temp(count, 2) = length(d);
    beginNend_temp(count, 3:4) = age_estimate(count,:);
    
    l(count,1) = length(d);
    
    if testM == true
%     t_axis = beginNend_temp(count,3) : (beginNend_temp(count,4)-beginNend_temp(count,3))/(length(d)-1) : beginNend_temp(count,4);
%     plot(t_axis, d(:,3))
%     pause
    
        input = d;
        
        tB = beginNend_temp(count, 3);
        tE = beginNend_temp(count, 4);
        input_scaled = (tB+(input(:,1)-input(1,1))*(tE-tB)/(input(end,1)-input(1,1)));
        input_scaled_l_temp(count) = length(input_scaled);
        sortedDiff = sort(diff(input_scaled));
        targetDistance_temp = max(sortedDiff(ceil(length(sortedDiff)*0.1))/4, 1/10);
        discrete_target = tB : targetDistance_temp : tE;
        discrete_target_l_temp(count) = length(discrete_target);
        
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
        count_log_f = 0;
        for ilog = 1 : size(log_f,1)
            for jlog = 1 : size(log_f,2)
%                 log_f{ilog,jlog} = zeros(4*dis_bin_input(ilog),1);
                count_log_f = count_log_f + 4*dis_bin_input(ilog);
            end
        end
%         whos_log_f = whos('log_f');
%         size_log_f_temp(count,1) = whos_log_f.bytes;        
        size_log_f_temp(count,1) = length(input_scaled)*length(discrete_target)*112 + count_log_f*8;
        
    end
end


[ll, li] = sort(l,1,'descend');
for i = 1 : length(ll)
    inputNAME{i} = inputNAME_temp{li(i),1};
    beginNend(i,:) = beginNend_temp(li(i),:);
    
    if testM == true
        discrete_target_l(i,1) = discrete_target_l_temp(li(i));
        input_scaled_l(i,1) = input_scaled_l_temp(li(i));
        
        maxbin(i,1) = maxbin_temp(li(i));
        maxbin_2(i,1) = maxbin_temp_2(li(i));
        size_log_f(i,1) = size_log_f_temp(li(i))/1024/1024/1024;
    end
end


if testM == true
    for i = 1 : length(ll)
        a(i,1) = max(maxbin_2(i), maxbin(i)*4)*discrete_target_l(i)*input_scaled_l(i)*8/1024/1024/1024;
        b(i,1) = max(maxbin(i)*4)*discrete_target_l(i)*input_scaled_l(i)*8/1024/1024/1024;
    end
    
    c = a - b;
    [cc ci] = sort(c, 'descend');
end
