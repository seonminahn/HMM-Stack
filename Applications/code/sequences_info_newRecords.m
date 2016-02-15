function [inputNAME, beginNend] = sequences_info_newRecords(dataL)
% clear all

% dataL = 1000;
testM = true;

% cd /Users/seonminahn/Dropbox/Brown_University/ALIGN/ProfileHMM/NewBenthicdata/NewRecords
oceanN = '../NewBenthicdata/NewRecords/';
pList = dir([oceanN '*.txt']);

beginD = 0;
endD = dataL;
count = 0;
precount = 0;
for i = 1 : length(pList)
    %     i
    %     disp([[oceanN pList(i).name]])
    d = load([oceanN pList(i).name]);
    
    clear seq
    [seq(:,1),label]=unique(d(:,1));
    seq(:,2)=d(label,2);
    
    begin1 = find(seq(:,2)>=beginD, 1, 'first');
    switch pList(i).name
        % Age unit is kyr
        case {'AII107-131.txt', 'MD03-2692.txt', 'MD03-2705.txt'}
            end1 = find(seq(:,2)<=endD*1000, 1, 'last');
        otherwise
            end1 = find(seq(:,2)<=endD, 1, 'last');
    end
    
    if isempty(end1)==0 || sum(isnan(seq(:,2)))>0
        switch pList(i).name
            case {'AII107-131.txt', 'MD03-2692.txt', 'MD03-2705.txt'}
                % Age unit is kyr
                count = count + 1;
                inputNAME_temp{count,1} = [oceanN, pList(i).name];
                
                beginNend_temp(count, 1) = begin1;              % Column index of the first point
                beginNend_temp(count, 2) = end1;                % Column index of the last point
                beginNend_temp(count, 3) = floor(seq(begin1,2)/1000);  % Age estimate of the first point
                beginNend_temp(count, 4) =  ceil(seq(end1,  2)/1000);     % Age estimate of the last point
                l(count,1) = beginNend_temp(count, 2) - beginNend_temp(count, 1) + 1;
                
            case {'C9001C.txt', 'U1313.txt', 'GIK23415.txt', 'ODP1147.txt'}
                % Skip
                % C9001C.txt, GIK23415.txt: No informtation
                % U1313.txt: Not sorted
                % ODP1147.txt: Age estimates from MetaData.xlsx: 132-32720? 3272?
                
            case 'DSDP594.txt'
                % Age estimates from MetaData.xlsx: 4-1000
                % some missing age estimates up to 162, then all NaNs in txt
                if endD >= 1000
                    count = count + 1;
                    inputNAME_temp{count,1} = [oceanN, pList(i).name];
                    
                    beginNend_temp(count, 1) = 1;
                    beginNend_temp(count, 2) = length(seq);
                    beginNend_temp(count, 3) = 4;
                    beginNend_temp(count, 4) = 1000;
                    l(count,1) = beginNend_temp(count, 2) - beginNend_temp(count, 1) + 1;
                end
                
            case 'GeoB3388.txt'
                % Age estimates from MetaData.xlsx: 10-1100
                % No age estimate in txt
                if endD >= 1100
                    count = count + 1;
                    inputNAME_temp{count,1} = [oceanN, pList(i).name];
                    
                    beginNend_temp(count, 1) = 1;
                    beginNend_temp(count, 2) = length(seq);
                    beginNend_temp(count, 3) = 10;
                    beginNend_temp(count, 4) = 1100;
                    l(count,1) = beginNend_temp(count, 2) - beginNend_temp(count, 1) + 1;
                end
                
            case 'GIK12392.txt'
                % Age estimates from MetaData.xlsx: 0-150
                % No age estimate in txt
                if endD >= 150
                    count = count + 1;
                    inputNAME_temp{count,1} = [oceanN, pList(i).name];
                    
                    beginNend_temp(count, 1) = 1;
                    beginNend_temp(count, 2) = length(seq);
                    beginNend_temp(count, 3) = 0;
                    beginNend_temp(count, 4) = 150;
                    l(count,1) = beginNend_temp(count, 2) - beginNend_temp(count, 1) + 1;
                end
                
            case 'GIK13519.txt'
                % Age estimates from MetaData.xlsx: 0-750
                % some missing age estimates up to 139, then all NaNs in txt
                if endD >= 750
                    count = count + 1;
                    inputNAME_temp{count,1} = [oceanN, pList(i).name];
                    
                    beginNend_temp(count, 1) = 1;
                    beginNend_temp(count, 2) = length(seq);
                    beginNend_temp(count, 3) = 0;
                    beginNend_temp(count, 4) = 750;
                    l(count,1) = beginNend_temp(count, 2) - beginNend_temp(count, 1) + 1;
                end
                
            case 'GIK23414.txt'
                % Age estimates from MetaData.xlsx: 0-250
                % No age estimate in txt
                if endD >= 250
                    count = count + 1;
                    inputNAME_temp{count,1} = [oceanN, pList(i).name];
                    
                    beginNend_temp(count, 1) = 1;
                    beginNend_temp(count, 2) = length(seq);
                    beginNend_temp(count, 3) = 0;
                    beginNend_temp(count, 4) = 250;
                    l(count,1) = beginNend_temp(count, 2) - beginNend_temp(count, 1) + 1;
                end
                
            case 'M35027.txt'
                % Age estimates from MetaData.xlsx: 10-350
                % No age estimate in txt
                if endD >= 350
                    count = count + 1;
                    inputNAME_temp{count,1} = [oceanN, pList(i).name];
                    
                    beginNend_temp(count, 1) = 1;
                    beginNend_temp(count, 2) = length(seq);
                    beginNend_temp(count, 3) = 10;
                    beginNend_temp(count, 4) = 350;
                    l(count,1) = beginNend_temp(count, 2) - beginNend_temp(count, 1) + 1;
                end
                
            case 'MD01-2416.txt'
                % Age estimates from MetaData.xlsx: 0-450
                % No age estimate in txt
                if endD >= 450
                    count = count + 1;
                    inputNAME_temp{count,1} = [oceanN, pList(i).name];
                    
                    beginNend_temp(count, 1) = 1;
                    beginNend_temp(count, 2) = length(seq);
                    beginNend_temp(count, 3) = 0;
                    beginNend_temp(count, 4) = 450;
                    l(count,1) = beginNend_temp(count, 2) - beginNend_temp(count, 1) + 1;
                end
                
            case 'MD96-2080.txt'
                % Age estimates from MetaData.xlsx: 0-450
                % No age estimate in txt
                if endD >= 450
                    count = count + 1;
                    inputNAME_temp{count,1} = [oceanN, pList(i).name];
                    
                    beginNend_temp(count, 1) = 1;
                    beginNend_temp(count, 2) = length(seq);
                    beginNend_temp(count, 3) = 0;
                    beginNend_temp(count, 4) = 450;
                    l(count,1) = beginNend_temp(count, 2) - beginNend_temp(count, 1) + 1;
                end
                
            case 'MD02-2575.txt'
                % Need to modify data
                % Duplicated data and diffrent unit in age
                % Age estimates from MetaData.xlsx: 0-400
                if endD >= 400
                    count = count + 1;
                    inputNAME_temp{count,1} = [oceanN, pList(i).name];
                    
                    beginNend_temp(count, 1) = 1;
                    beginNend_temp(count, 2) = length(seq);
                    beginNend_temp(count, 3) = 0;
                    beginNend_temp(count, 4) = 400;
                    l(count,1) = beginNend_temp(count, 2) - beginNend_temp(count, 1) + 1;
                end
                
            case 'ODP723.txt'
                % Age estimates from MetaData.xlsx: 0-1500
                % No age estimate in txt
                if endD >= 1500
                    count = count + 1;
                    inputNAME_temp{count,1} = [oceanN, pList(i).name];
                    
                    beginNend_temp(count, 1) = 1;
                    beginNend_temp(count, 2) = length(seq);
                    beginNend_temp(count, 3) = 0;
                    beginNend_temp(count, 4) = 1500;
                    l(count,1) = beginNend_temp(count, 2) - beginNend_temp(count, 1) + 1;
                end
                
            case 'ODP1209.txt'
                % Age estimates from MetaData.xlsx: 0-350
                % No age estimate in txt
                if endD >= 350
                    count = count + 1;
                    inputNAME_temp{count,1} = [oceanN, pList(i).name];
                    
                    beginNend_temp(count, 1) = 1;
                    beginNend_temp(count, 2) = length(seq);
                    beginNend_temp(count, 3) = 0;
                    beginNend_temp(count, 4) = 350;
                    l(count,1) = beginNend_temp(count, 2) - beginNend_temp(count, 1) + 1;
                end
                
            case 'ODP1236.txt'
                % Age estimates from MetaData.xlsx: 2500-5300
                % No age estimate in txt
                if endD >= 5300
                    count = count + 1;
                    inputNAME_temp{count,1} = [oceanN, pList(i).name];
                    
                    beginNend_temp(count, 1) = 1;
                    beginNend_temp(count, 2) = length(seq);
                    beginNend_temp(count, 3) = 2500;
                    beginNend_temp(count, 4) = 5300;
                    l(count,1) = beginNend_temp(count, 2) - beginNend_temp(count, 1) + 1;
                end
                
            case 'ODP1237.txt'
                % Age estimates from MetaData.xlsx: 4000-5000
                % No age estimate in txt
                if endD >= 5000
                    count = count + 1;
                    inputNAME_temp{count,1} = [oceanN, pList(i).name];
                    
                    beginNend_temp(count, 1) = 1;
                    beginNend_temp(count, 2) = length(seq);
                    beginNend_temp(count, 3) = 4000;
                    beginNend_temp(count, 4) = 5000;
                    l(count,1) = beginNend_temp(count, 2) - beginNend_temp(count, 1) + 1;
                end
                
            case 'ODP1239.txt'
                % Age estimates from MetaData.xlsx: 2700-4900
                % No age estimate in txt
                if endD >= 4900
                    count = count + 1;
                    inputNAME_temp{count,1} = [oceanN, pList(i).name];
                    
                    beginNend_temp(count, 1) = 1;
                    beginNend_temp(count, 2) = length(seq);
                    beginNend_temp(count, 3) = 2700;
                    beginNend_temp(count, 4) = 4900;
                    l(count,1) = beginNend_temp(count, 2) - beginNend_temp(count, 1) + 1;
                end
                
            case 'ODP1241.txt'
                % Age estimates from MetaData.xlsx: 0-5691
                % some missing age estimates up to 2129, then all NaNs in txt
                if endD <= 5690
                    count = count + 1;
                    inputNAME_temp{count,1} = [oceanN, pList(i).name];
                    
                    beginNend_temp(count, 1) = begin1;
                    beginNend_temp(count, 2) = end1;
                    beginNend_temp(count, 3) = floor(seq(begin1,2));
                    beginNend_temp(count, 4) =  ceil(seq(end1,2));
                    l(count,1) = beginNend_temp(count, 2) - beginNend_temp(count, 1) + 1;
                    
                elseif endD > 5690
                    count = count + 1;
                    inputNAME_temp{count,1} = [oceanN, pList(i).name];
                    
                    beginNend_temp(count, 1) = 1;
                    beginNend_temp(count, 2) = length(seq);
                    beginNend_temp(count, 3) = 0;
                    beginNend_temp(count, 4) = 5690;
                    l(count,1) = beginNend_temp(count, 2) - beginNend_temp(count, 1) + 1;
                end
                
                
            case 'RC13-229.txt'
                % Age estimates from MetaData.xlsx: 0-750
                % some missing age estimates up to 32, then all NaNs in txt
                if endD >= 750
                    count = count + 1;
                    inputNAME_temp{count,1} = [oceanN, pList(i).name];
                    
                    beginNend_temp(count, 1) = 1;
                    beginNend_temp(count, 2) = length(seq);
                    beginNend_temp(count, 3) = 0;
                    beginNend_temp(count, 4) = 750;
                    l(count,1) = beginNend_temp(count, 2) - beginNend_temp(count, 1) + 1;
                end
                
            case 'V29-135.txt'
                % Age estimates from MetaData.xlsx: 0-300
                % some missing age estimates up to 58, then all NaNs in txt
                if endD >= 300
                    count = count + 1;
                    inputNAME_temp{count,1} = [oceanN, pList(i).name];
                    
                    beginNend_temp(count, 1) = 1;
                    beginNend_temp(count, 2) = length(seq);
                    beginNend_temp(count, 3) = 0;
                    beginNend_temp(count, 4) = 300;
                    l(count,1) = beginNend_temp(count, 2) - beginNend_temp(count, 1) + 1;
                end
                
            case 'V29-202.txt'
                % Age estimates from MetaData.xlsx: 0-200
                % No age estimate in txt
                if endD >= 200
                    count = count + 1;
                    inputNAME_temp{count,1} = [oceanN, pList(i).name];
                    
                    beginNend_temp(count, 1) = 1;
                    beginNend_temp(count, 2) = length(seq);
                    beginNend_temp(count, 3) = 0;
                    beginNend_temp(count, 4) = 200;
                    l(count,1) = beginNend_temp(count, 2) - beginNend_temp(count, 1) + 1;
                end
                
            case 'V30-97.txt'
                % Age estimates from MetaData.xlsx: 0-250
                % No age estimate in txt
                if endD >= 250
                    count = count + 1;
                    inputNAME_temp{count,1} = [oceanN, pList(i).name];
                    l(count,1) = length(seq);
                    
                    beginNend_temp(count, 1) = 1;
                    beginNend_temp(count, 2) = length(seq);
                    beginNend_temp(count, 3) = 0;
                    beginNend_temp(count, 4) = 250;
                    l(count,1) = beginNend_temp(count, 2) - beginNend_temp(count, 1) + 1;
                end
                
            otherwise
                % All other records
                % 'MD12-3409.txt': Negative age estimate at the first point
                % 'MV0502.txt': More data after 2975, but not used
                count = count + 1
                inputNAME_temp{count,1} = [oceanN, pList(i).name];
                
                beginNend_temp(count, 1) = begin1;              % Column index of the first point
                beginNend_temp(count, 2) = end1;                % Column index of the last point
                beginNend_temp(count, 3) = floor(seq(begin1,2));  % Age estimate of the first point
                beginNend_temp(count, 4) =  ceil(seq(end1,2));     % Age estimate of the last point
                l(count,1) = beginNend_temp(count, 2) - beginNend_temp(count, 1) + 1;
                
        end
    end
    
    if testM == true
        if precount ~= count
            input = seq(beginNend_temp(count, 1):beginNend_temp(count, 2),:);
            
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
        end
        precount = count;
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
