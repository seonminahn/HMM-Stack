%% forward_point_fine;
% log_f=zeros(length(input_scaled_new),length(discrete_target),length(discrete_target));
%%%SA%%%
log_f = cell(length(input_scaled_new),length(discrete_target));
for i = 1 : size(log_f,1)
    for j = 1 : size(log_f,2)
        log_f{i,j} = zeros(4*dis_bin_input_new(i),1);
    end
end
% log_f = zeros(length(input_scaled_new),length(discrete_target), prac_max_dis_bin_input_new);
% log_f = zeros(length(input_scaled_new),length(discrete_target), 4*max_dis_bin_input_new);
%%%%%%%%

%%%SA%%% mean_tar should be defined outside spmd for profile HMM,
%%%computed inside sequnces_info 1
% mean_tar=mean(target(:,2));
%%%%%%%%

%initiliazation

% for i=1:length(discrete_target)
%     %     log_f(:,1:i,i) = -inf;
%     for j=1:i-1
%         log_f(:,j,i)=-inf;
%     end
% end

%setting permitable range of mapping rates
% mm=4;
%
% for i=2:length(input_scaled_new)
%     for j = 1:length(discrete_target)
%         if j+(dis_bin_input_new(i))*mm+1<=length(discrete_target)
%             log_f(i,min(j+(dis_bin_input_new(i))*mm+1,length(discrete_target)):end,j) = -inf;
%         end
%     end
% end

for i=2:size(log_f,1)
    for j=1:size(log_f,2)
        if j<=4*dis_bin_input_new(i)
            for k=1:j
                
                rr = (dis_bin_input_new(i)-1)/max(k-1,0.5);
                
                if rr<0.25 || rr>4
                    log_f{i,j}(k)=-inf;
                end
            end
            for k=j+1:4*dis_bin_input_new(i)
                log_f{i,j}(k)=-inf;
            end
        else
            for k=1:4*dis_bin_input_new(i)
                rr = (dis_bin_input_new(i)-1)/max(k-1,0.5);
                if rr<0.25 || rr>4
                    log_f{i,j}(k)=-inf;
                end
            end
            
        end
    end
end

% for j=1:size(log_f,2)
%     for k=j+1:size(log_f,3)
%         log_f(1,j,k)=-inf;
%     end
% end


log_f{1,1}(1:end) = log_pi_begin(2)+log_grid(1,1);

%less than gan_con_tar points gap in target

for i = 2:gap_con_tar
    log_f{1,i}(1:end) = log_grid(1,i)+log_pi_begin(3)+log_tao_y_begin(1)*(i-2)+ log_tao_y_begin(2);
end

for i=gap_con_tar+1:length(discrete_target)
    log_f{1,i}(1:end)=-inf;
end

for i=2:gap_con
    log_f{i,1}(:) = log_grid(i,1)+log_pi_begin(1)+log_tao_x_begin(1)*(i-2)+log_tao_x_begin(2);
    for j=1:i-1
        %%%SA%%%
        % log_f(i,1,:) = log_f(i,1,:)+ log(normpdf(input(j,2)-mean_tar,mu,sigma));
        log_f{i,1}(:) = log_f{i,1}(:)+ log(normpdf(input(j,2)-mean_tar, mean_mu, mean_sigma));
        %%%%%%%%
    end
end

for i=gap_con+1:length(input_scaled_new)
    log_f{i,1}(:)=-inf;
end


%setup for the second point if first point is mapped
for i=2:length(discrete_target)
    for j=1:4*dis_bin_input_new(2)
        if log_f{2,i}(j)==0
            
            rr = (dis_bin_input_new(2)-1)/max(j-1,0.5);
            
            if rr<min_r || rr>max_r
                log_f{2,i}(j)=-inf;
            else
                nearest_index = ratio_track(dis_bin_input_new(2),j);
                log_f{2,i}(j) = log_f{1,i-j+1}(1)+log_grid(2,i)...
                    +log_delta_match(1)...   %SA% log_delta_r(1)
                    +log(pi_0(nearest_index));
            end
        end
        
    end
end


for i=2:length(input_scaled_new)
    %%%SA%%%%
    for j=2:min(4*dis_bin_input_new(i), size(log_f,2))
%     for j=2:4*dis_bin_input_new(i)
    %%%%%%%%
        if log_f{i,j}(j)==0
            
            rr = ratio_track_new(i,j);
            
            if rr<min_r || rr>max_r
                log_f{i,j}(j)=-inf;
            else
                nearest_index = ratio_track(dis_bin_input_new(i),j);
                log_f{i,j}(j) = log_f{i-1,1}(1)+log_grid(i,j)...
                    +log_delta_match(1)...   %SA% log_delta_r(1)
                    + log(pi_0(nearest_index));
            end
        end
    end
end


%recursion
% mem_use()/1024
ssum = 0;
for i=3:length(input_scaled_new)
    for j=1:length(discrete_target)
        j2_range = find(log_f{i,j}(:) == 0);
        for j2 = j2_range'
%         for j2=1:4*dis_bin_input_new(i)
%             if log_f{i,j}(j2)==0
                rr = ratio_track_new(i,j2);
                
                log_f{i,j}(j2) = log_grid(i,j);
                [max_logf,max_ind]=max(log_f{i-1,j-j2+1}(:));
                if max_logf==-inf
                    log_f{i,j}(j2) = -inf;
                else
                    pre_rr = ratio_track_new(i-1,max_ind);
                    
                    if pre_rr<0.9220
                        qq1=1;
                    elseif pre_rr>1.085
                        qq1=3;
                    else
                        qq1=2;
                    end
                    
                    if rr<0.9220
                        qq2=1;
                    elseif rr>1.085
                        qq2=3;
                    else
                        qq2=2;
                    end
                    
                    ssum = group_transition(qq1,qq2);

                    for t=1:4*dis_bin_input_new(i-1)
                        if log_f{i-1,j-j2+1}(t)~=-inf && t~=max_ind
                            
                            cc = ratio_track_new(i-1,t);
                            
                            if cc<0.9220
                                qq3=1;
                            elseif cc>1.085
                                qq3=3;
                            else
                                qq3=2;
                            end
                            ssum = ssum+group_transition(qq3,qq2)*exp(log_f{i-1,j-j2+1}(t)-max_logf);
                            
                        end;
                    end;
                    
%                     t_range = find(log_f{i-1,j-j2+1}(:) ~= -inf);
%                     t_range(t_range == max_ind) = [];
%                     if isempty(t_range)
%                     else                        
%                         cc = ratio_track_new(i-1,t_range);
%                         
%                         qq3 = 2*ones(size(cc));
%                         qq3(cc<0.9220) = 1;
%                         qq3(cc>1.085) = 3;
%                         
%                         ssum = ssum + group_transition(qq3,qq2)'*exp(log_f{i-1,j-j2+1}(t_range)-max_logf);
%                     end;
                    
                    

                                                
                    ssum = ssum*transition_prob_track(i,j2);
                    %%%SA%%%
                    % log_f{i,j}(j2)=log_f{i,j}(j2)+max_logf+log(ssum)+log_delta_r(1);
                    log_f{i,j}(j2)=log_f{i,j}(j2)+max_logf+log(ssum)+log_delta_match(1);
                    %%%%%%%%
                end
%             end
        end
    end
end


% tic
% for i=3:length(input_scaled_new)
%   % i
%
%     for j=1:length(discrete_target)
%        % j
%         for j2=max(2,j-4*dis_bin_input_new(i)):j
%             if log_f(i,j,j2)==0
%
%                     rr = ratio_track_new(i,j-j2+1);
%
%
%                     log_f(i,j,j2) = log_grid(i,j);
%                     [max_logf,max_ind]=max(log_f(i-1,j2,1:j2));
%                     if max_logf==-inf
%                         log_f(i,j,j2) = -inf;
%                     else
%
%                        pre_rr = ratio_track_new(i-1,j2-max_ind+1);
%
%                         if pre_rr<0.9220
%                             qq1=1;
%                         elseif pre_rr>1.085
%                             qq1=3;
%                         else
%                             qq1=2;
%                         end
%
%                         if rr<0.9220
%                             qq2=1;
%                         elseif rr>1.085
%                             qq2=3;
%                         else
%                             qq2=2;
%                         end
%
%                         ssum = group_transition(qq1,qq2);
%                         for t=1:j2
%                             if log_f(i-1,j2,t)~=-inf && t~=max_ind
%
%                                 cc = ratio_track_new(i-1,j2-t+1);
%
%                                 if cc<0.9220
%                                     qq3=1;
%                                 elseif cc>1.085
%                                     qq3=3;
%                                 else
%                                     qq3=2;
%                                 end
%
%
%
%
%                                 ssum = ssum+group_transition(qq3,qq2)*exp(log_f(i-1,j2,t)-max_logf);
%
%                             end;
%                         end;
%                         ssum = ssum*transition_prob_track(i,j-j2+1);
%                         log_f(i,j,j2)=log_f(i,j,j2)+max_logf+log(ssum)+log_delta_r(1);
%                     end
%
%             end
%         end
%
%     end
%
% end
% toc


%add nomatch penalty for ending parts
%existing unmatched bin in target
gap_add_f_x = zeros(length(input_scaled_new)-1,length(discrete_target)); %f(i,N2,k)
gap_add_f_y = zeros(length(discrete_target)-1,length(discrete_target)-1); %f(N1,j,k)
gap_add_f_e = zeros(length(discrete_target),1);

for i=1:length(input_scaled_new)-1
    if i<=length(input_scaled_new)-gap_con
        gap_add_f_x(i,:)=-inf;
    else
        for k=1:length(discrete_target)
            if length(discrete_target)-k+1>0 && length(discrete_target)-k+1<=4*dis_bin_input_new(i)
                if log_f{i,length(discrete_target)}(length(discrete_target)-k+1)==-inf
                    gap_add_f_x(i,k)=-inf;
                else
    
% %SA%                    gap_add_f_x(i,k) = log_f{i,length(discrete_target)}(length(discrete_target)-k+1)+log_delta_r(2)+log_tao_x_end(1)*...
                    gap_add_f_x(i,k) = log_f{i,length(discrete_target)}(length(discrete_target)-k+1)+log_delta_match(2)+log_delta_end(2)+log_tao_x_end(1)*...
                        (length(input_scaled_new)-i-1)+log_tao_x_end(2);
                    for j=i+1:length(input_scaled_new)
                        %%%SA%%%
                        gap_add_f_x(i,k)= gap_add_f_x(i,k)+log(normpdf(input_scaled_new(j,2)-mean_tar,mean_mu,mean_sigma));
                        % gap_add_f_x(i,k)= gap_add_f_x(i,k)+log(normpdf(input_scaled_new(j,2)-mean_tar,mu,sigma));
                        %%%%%%%%
                    end
                end
            else
                gap_add_f_x(i,k)=-inf;
            end
        end
    end
end

for j=1:length(discrete_target)-1
    if j<length(discrete_target)-gap_con_tar
        gap_add_f_y(j,:)=-inf;
    else
        for k=1:length(discrete_target)-1
            if j-k+1>4*dis_bin_input_new(end) || j-k+1<1
                gap_add_f_y(j,k)=-inf;
            elseif log_f{length(input_scaled_new),j}(j-k+1)==-inf
                gap_add_f_y(j,k)=-inf;
            else
% %SA%             gap_add_f_y(j,k) = log_f{length(input_scaled_new),j}(j-k+1)+log_delta_r(3)+log_tao_y_end(1)*...
                gap_add_f_y(j,k) = log_f{length(input_scaled_new),j}(j-k+1)+log_delta_match(2)+log_delta_end(3)+log_tao_y_end(1)*...
                    (length(discrete_target)-j-1)+log_tao_y_end(2);
            end
        end
    end
end

for i=1:length(discrete_target)
    if length(discrete_target)-i+1>0 && length(discrete_target)-i+1<=4*dis_bin_input_new(end)
% %SA%        gap_add_f_e(i) = log_f{length(input_scaled_new),length(discrete_target)}(length(discrete_target)-i+1)+log_delta_r(3);
        gap_add_f_e(i) = log_f{length(input_scaled_new),length(discrete_target)}(length(discrete_target)-i+1)+log_delta_match(2)+log_delta_end(1);
    else
        gap_add_f_e(i)=-inf;
    end
end


%termination
%log_fE is just the log likelihood of the emission sequence
tt1 = max(gap_add_f_e);
tt2 = max(max(gap_add_f_x));
tt3 = max(max(gap_add_f_y));
tt = max([tt1 tt2 tt3]);

log_fE=0;
for s = 1:length(discrete_target)
    if gap_add_f_e(s)~=-inf
        log_fE = log_fE + exp(gap_add_f_e(s)-tt);
    end
end
for i=1:length(input_scaled_new)-1
    for k=1:length(discrete_target)
        if gap_add_f_x(i,k)~=-inf
            log_fE = log_fE + exp(gap_add_f_x(i,k)-tt);
        end
    end
end
for i=1:length(discrete_target)-1
    for j=1:length(discrete_target)-1
        if gap_add_f_y(i,j)~=-inf
            log_fE = log_fE + exp(gap_add_f_y(i,j)-tt);
        end
    end
end

log_fE=tt+log(log_fE);

fprintf('Forward algorithm is done.\t\t')
toc