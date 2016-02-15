%%%SA%%%
% load(['../../ProfileHMM_results/26-May-2015/LR04plusNew_6000_IC1_step_norm0_input', num2str(selectInput), '_bb7.mat'])
% load(['../../ProfileHMM_results/16-Sep-2015/LR04_4000_IC1_step_norm0_input', num2str(selectInput), '_bb4.mat'])
% resultNS = [ coreNAME '_' num2str(dataL) '_IC' num2str(exactStart) '_' est_tag '_norm' num2str(dataNormalize) '_input', num2str(selectInput), '_bb'];
% resultS = ['../../ProfileHMM_results/' r_date '/' resultNS '*' ];
dList = dir(resultS);
for i = 1 : length(dList)
    count(i) = sscanf(dList(i).name, [resultNS '%d']);
end
loadFN = [strrep(resultS, '*', num2str(max(count)))];
disp(' ')
disp(['Loading ' loadFN])
load([loadFN '.mat'])
disp('Finding median age estimates...')
calTIC = tic;
%%%%%%%%

%%
%given samples, return centroid alignment
% - cen_seq: centroid alignment with aligned bin index info
% - cen_ali: centroid alignment with aligned age info
MarginalPairCount = zeros(size(sample_ali,1),length(discrete_target)+2);

for s = 1:size(sample_seq,2)
    vv=find(sample_seq(:,s)~=0);
    if vv(1)~=1
        for i=1:vv(1)-1
            MarginalPairCount(i,length(discrete_target)+1)= MarginalPairCount(i...
                ,length(discrete_target)+1)+1;
        end
    end
    if vv(end)~=size(sample_seq,1)
        for i=vv(end)+1:size(sample_seq,1)
            MarginalPairCount(i,length(discrete_target)+2)=MarginalPairCount(i,...
                length(discrete_target)+2)+1;
        end
    end
    for t1=vv(1):vv(end)
        MarginalPairCount(t1,sample_seq(t1,s))=MarginalPairCount(t1,sample_seq(t1,s))+1;
        
    end
end

ProbMarginalPair =zeros(size(sample_ali,1),length(discrete_target)+2);
for i=1:size(ProbMarginalPair,1)
    if sum(MarginalPairCount(i,:))==0
        ProbMarginalPair(i,:)=0;
    else
        ProbMarginalPair(i,:) = MarginalPairCount(i,:)/sum(MarginalPairCount(i,:));
    end
end


ExpAbsDiff = zeros(size(sample_ali,1),length(discrete_target)+2);

for i=1:size(sample_ali,1)
    for j=1:length(discrete_target)+2
        if j<=length(discrete_target)
            s=0;
            for tt=1:length(discrete_target)
                s=s+abs(tt-j)*ProbMarginalPair(i,tt);
            end
            ExpAbsDiff(i,j)=s;
        elseif j==length(discrete_target)
            s=0;
            for tt=1:length(discrete_target)
                s=s+abs(tt-0)*ProbMarginalPair(i,tt);
            end
            ExpAbsDiff(i,j)=s;
        else
            s=0;
            for tt=1:length(discrete_target)
                s=s+abs(tt-length(discrete_target))*ProbMarginalPair(i,tt);
            end
            ExpAbsDiff(i,j)=s;
        end
        
    end
end

ExpAbsDiff=-ExpAbsDiff;

multip=4;
%%%SA%%%
% VV = zeros(size(sample_seq,1),length(discrete_target),4*max(dis_bin_input_new));
VV = cell(size(sample_seq,1),length(discrete_target));
for i = 1 : size(VV,1)
    for j = 1 : size(VV,2)
        VV{i,j} = zeros(4*dis_bin_input_new(i),1);
    end
end
% All lines with VV are modified according the new definition
%%%%%%%%

for i=2:size(VV,1)
    for j=1:size(VV,2)
        if j<=size(VV{i,j},1)
            for k=1:j
                rr = (dis_bin_input_new(i)-1)/max(k-1,0.5);
                if rr<0.25 || rr>4
                    VV{i,j}(k)=-inf;
                end
            end
            for k=j+1:size(VV{i,j},1)
                VV{i,j}(k)=-inf;
            end
        else
            for k=1:size(VV{i,j},1)
                rr = (dis_bin_input_new(i)-1)/max(k-1,0.5);
                if rr<0.25 || rr>4
                    VV{i,j}(k)=-inf;
                end
            end
            
        end
    end
end



% for i=1:length(discrete_target)
%     for j=1:i-1
%         VV(:,j,i)=-inf;
%     end
% end


% for i=2:size(sample_seq,1)
%     for j = i:length(discrete_target)
%         if j+(dis_bin_input_new(i))*multip+1<=length(discrete_target)
%             VV(i,j+(dis_bin_input_new(i))*multip+1:end,j) = -inf;
%         end
%     end
% end

for i=2:size(VV,1)
    % i
    for j=1:size(VV,2)
        if j<=size(VV{i,j},1)
            for k=1:j
                %  rr = ratio_track_new(i,j-k+1);
                
                rr = (dis_bin_input_new(i)-1)/max(k-1,0.5);
                if rr<0.25 || rr>4
                    VV{i,j}(k)=-inf;
                end
            end
            for k=j+1:size(VV{i,j},1)
                VV{i,j}(k)=-inf;
            end
        else
            for k=1:size(VV{i,j},1)
                rr=(dis_bin_input_new(i)-1)/max(k-1,0.5);
                if rr<0.25 || rr>4
                    VV{i,j}(k)=-inf;
                end
            end
            
        end
    end
end

%log_v(1,1,1) = log(normpdf(input_scaled(1,2)-bin_mean(1),0,sigma));
VV{1,1}(:) = ExpAbsDiff(1,1);
for i = 2:gap_con_tar
    VV{1,i}(:) = ExpAbsDiff(1,i);
end
for i=gap_con_tar+1:length(discrete_target)
    VV{1,i}(:)=-inf;
end

for i=2:gap_con
    VV{i,1}(:) = sum(ExpAbsDiff(1:i-1,length(discrete_target)+1))...
        +ExpAbsDiff(i,1);
end
for i=gap_con+1:length(input_scaled_new)
    VV{i,1}(:)=-inf;
end


%recursion

for i=2:length(input_scaled_new)
    %  i
    for j=1:length(discrete_target)
        for j2=1:j
            tmp=j-j2+1;
            if tmp>0 && tmp<=size(VV{i,j},1)
                if VV{i,j}(tmp)==0
                    
                    %  if log_v(i,j,j2)==-inf
                    %   VV(i,j,j2)=-inf;
                    % else
                    VV{i,j}(tmp) =ExpAbsDiff(i,j)+max(VV{i-1,j2}(:));
                    
                    
                    % end
                end
            end
        end
    end
end

gap_add_v_x = zeros(length(input_scaled_new)-1,length(discrete_target)); %f(i,N2,k)
gap_add_v_y = zeros(length(discrete_target)-1,length(discrete_target)-1); %f(N1,j,k)
gap_add_v_e = zeros(length(discrete_target),1);


for i=1:length(input_scaled_new)-1
    if i<length(input_scaled_new)-gap_con
        gap_add_v_x(i,:)=-inf;
    else
        for k=1:length(discrete_target)
            if length(discrete_target)-k+1>0 && length(discrete_target)-k+1<=size(VV{i,1},1)
                if VV{i,length(discrete_target)}(length(discrete_target)-k+1)==-inf
                    gap_add_v_x(i,k)=-inf;
                else
                    gap_add_v_x(i,k) = VV{i,length(discrete_target)}(length(discrete_target)-k+1)+sum(ExpAbsDiff(i+1:...
                        length(input_scaled_new),end));
                    
                end
            else
                gap_add_v_x(i,k)=-inf;
            end
            
        end
    end
end
for j=1:length(discrete_target)-1
    if j<length(discrete_target)-gap_con_tar
        gap_add_v_y(j,:)=-inf;
    else
        for k=1:length(discrete_target)-1
            tmp=j-k+1;
            if tmp>0 && tmp<=size(VV{end,1},1)
                if VV{length(input_scaled_new),j}(tmp)==-inf
                    gap_add_v_y(j,k)=-inf;
                else
                    gap_add_v_y(j,k) = VV{length(input_scaled_new),j}(tmp);
                    
                end
            else
                gap_add_v_y(j,k)=-inf;
            end
        end
    end
end
for i=1:length(discrete_target)
    tmp = length(discrete_target)-i+1;
    if tmp>0 && tmp<=size(VV{end,1},1)
        gap_add_v_e(i) = VV{length(input_scaled_new),length(discrete_target)}(tmp);
    else
        gap_add_v_e(i)=-inf;
    end
end

%backtrace

ttt=-inf;
for i=1:length(input_scaled_new)-1
    for j=1:length(discrete_target)
        if gap_add_v_x(i,j)>ttt
            ttt=gap_add_v_x(i,j);
            ind_i=i;
            ind_j=length(discrete_target);
            ind_k=j;
        end
    end
end
for i=1:length(discrete_target)-1
    for j=1:length(discrete_target)-1
        if gap_add_v_y(i,j)>ttt
            ttt=gap_add_v_y(i,j);
            ind_i=length(input_scaled_new);
            ind_j=i;
            ind_k=j;
        end
    end
end
for i=1:length(discrete_target)
    if gap_add_v_e(i)>ttt
        ttt=gap_add_v_e(i);
        ind_i=length(input_scaled_new);
        ind_j=length(discrete_target);
        ind_k=i;
    end
end

%[ind_i ind_j ind_k] = optimal_ali_ending(log_v_gap);
median_seq = zeros(length(input_scaled_new),1);
median_seq(ind_i) = ind_j;
median_seq(ind_i-1) = ind_k;

for i = ind_i-2:-1:1
    %sss = zeros(ind_k,1);
    sss=zeros(ind_k+1-max(1,ind_k+1-size(VV{i,1},1)),1);
    
    for j = max(1,ind_k+1-size(VV{i+1,ind_k},1)):ind_k
        % nearest_index = ratio_track(dis_bin_input(i+2),cen_seq(i+2)-ind_k+1);
        % if log_v(i+1,ind_k,j)==-inf
        %    sss(j)=inf;
        % else
        
        
        sss(j-max(1,ind_k+1-size(VV{i+1,ind_k},1))+1)=abs(ExpAbsDiff(i+2,median_seq(i+2))+VV{i+1,ind_k}(ind_k-j+1)-...
            VV{i+2,median_seq(i+2)}(median_seq(i+2)-ind_k+1));
        
        %end
    end
    %sss(ind_k+1)=abs(ProbMarginalPair(i+2,cen_seq(i+2))+...
    % +ProbMarginalPair(i+1,length(discrete_target)+1)-VV(i+2,cen_seq(i+2),ind_k));
    [aa bb]= min(sss);
    if bb==ind_k+1
        median_seq(i)=0;
        break
    else
        median_seq(i)=bb+max(1,ind_k+1-size(VV{i+1,ind_k},1))-1;
        ind_k=bb+max(1,ind_k+1-size(VV{i+1,ind_k},1))-1;
        if bb+max(1,ind_k+1-size(VV{i+1,ind_k},1))-1==1
            break;
        end
    end
end

median_ali=zeros(length(input_scaled_new),1);
for i=1:length(input_scaled_new)
    if median_seq(i)~=0
        median_ali(i) = discrete_target(median_seq(i));
    end
end

aa_end = max(find(median_ali~=0));
if aa_end~=length(median_ali)
    for i=aa_end+1:length(median_ali)
        median_ali(i)=median_ali(aa_end);
    end
end

%%
%%%SA%%%
calT = toc(calTIC);
% save(['../../ProfileHMM_results/16-Sep-2015/LR04_4000_IC1_step_norm0_input', num2str(selectInput), '_bb4_median.mat'], 'median_seq', 'median_ali', 'calT')
save([loadFN '_median.mat'])
%%%%%%%%