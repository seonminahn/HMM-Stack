%% back_sample_point_fine;
%generate sample alignments -- back sampling
% - sample_seq contains the 1000 samples, each column represents one
% sample,similar for sample_ali
% - sample_seq: records the aligned bin index for alignments
% - sample_ale: records the aligned age info for alignments
sample_ali=zeros(length(input_scaled_new),1000);
sample_seq=zeros(length(input_scaled_new),1000);
%exp_logf=exp(log_f);



%% back_sampling_end_global;  (in back_sample_point_fine)
%sampling the end status
n_ali = length(discrete_target)+(length(input_scaled_new)-1)+(length(discrete_target)-1);
log_pp_ali = zeros(n_ali,1);
index_set = zeros(n_ali,2);

vv_ali=0;
vv2=0;
tt = find(gap_add_f_e~=-inf & isnan(gap_add_f_e)==0);
tt2=gap_add_f_e(tt);
tt3=max(tt2);
for i=1:length(tt)
    vv2=vv2+exp(tt2(i)-tt3);
end
vvv=log(vv2)+tt3;

% for i=1:length(discrete_target)
%     if gap_add_f_e(i)~=-inf && isnan(gap_add_f_e(i))==0
%         %vv_ali=vv_ali+1;
%         vv2 =vv2+ exp(gap_add_f_e(i));
%
%     end
% end

vv_ali=vv_ali+1;
index_set(vv_ali,1:2)=[length(input_scaled_new) length(discrete_target)];
log_pp_ali(vv_ali)=vvv;

sum_gap_add_f_x = zeros(length(input_scaled_new)-1,1);
for i=1:length(input_scaled_new)-1
    tt = max(gap_add_f_x(i,:));
    for j=1:length(discrete_target)
        sum_gap_add_f_x(i) = sum_gap_add_f_x(i) + exp(gap_add_f_x(i,j)-tt);
    end
    sum_gap_add_f_x(i) = log(sum_gap_add_f_x(i))+tt;
    if sum_gap_add_f_x(i)~=-inf && isnan(sum_gap_add_f_x(i))==0
        vv_ali = vv_ali+1;
        log_pp_ali(vv_ali) = sum_gap_add_f_x(i);
        index_set(vv_ali,:)=[i length(discrete_target)];
    end
end

sum_gap_add_f_y = zeros(length(discrete_target)-1,1);
for i=1:length(discrete_target)-1
    tt = max(gap_add_f_y(i,:));
    for j=1:length(discrete_target)-1
        sum_gap_add_f_y(i) = sum_gap_add_f_y(i) + exp(gap_add_f_y(i,j)-tt);
    end
    sum_gap_add_f_y(i) = log(sum_gap_add_f_y(i))+tt;
    if sum_gap_add_f_y(i)~=-inf && isnan(sum_gap_add_f_y(i))==0
        vv_ali = vv_ali+1;
        log_pp_ali(vv_ali) = sum_gap_add_f_y(i);
        index_set(vv_ali,:)=[length(input_scaled_new) i];
    end
end

log_pp_ali = log_pp_ali(1:vv_ali,:);
index_set = index_set(1:vv_ali,:);
end_set_sampling = zeros(1000,3);

for ss=1:1000
    cc = catigorical_random_log(log_pp_ali);
    
    ind_i=index_set(cc,1);
    ind_j=index_set(cc,2);
    
    if ind_i==length(input_scaled_new) && ind_j==length(discrete_target)
        cc2=catigorical_random_log(gap_add_f_e);
        ind_k=cc2;
    elseif ind_i==length(input_scaled_new)
        cc2 = catigorical_random_log(gap_add_f_y(ind_j,:));
        ind_k=cc2;
    else
        cc2 = catigorical_random_log(gap_add_f_x(ind_i,:));
        ind_k=cc2;
    end
    end_set_sampling(ss,:)=[ind_i ind_j ind_k];
end



%% BACK TO back_sample_point_fine;
for nn=1:size(sample_ali,2)
    curr_i=end_set_sampling(nn,1);
    curr_j=end_set_sampling(nn,2);
    s = end_set_sampling(nn,3);
    sample_seq(curr_i,nn) = curr_j;
    sample_seq(curr_i-1,nn) = s;
    %     tmp = curr_j-s+1;
    %     if dis_bin_input(curr_i)~=0 && tmp~=0
    %     pre_index = ratio_track(dis_bin_input(curr_i),tmp);
    %     elseif tmp~=0
    %         pre_index = ratio_track_0(tmp);
    %
    %     elseif dis_bin_input(curr_i)~=0
    %       pre_index = ratio_track_1(dis_bin_input(curr_i));
    %     else
    %         pre_index = 9;
    %     end
    pre_rr = ratio_track_new(curr_i, curr_j-s+1);
    %     if curr_j~=s
    %         pre_rr = dis_bin_input_new(curr_i)/(curr_j-s);
    %     else
    %         pre_rr=dis_bin_input_new(curr_i)/0.5;
    %     end
    curr_i=curr_i-1;
    curr_j=s;
    
    while curr_i~=1 && curr_j~=1
        u=rand;
        ind_nonzero = find(log_f{curr_i,curr_j}(:)~= -Inf);
        
        ssum2=zeros(length(ind_nonzero),1);
        pp=zeros(length(ind_nonzero),1);
        for tt=1:length(ind_nonzero)
            tmp = curr_j-ind_nonzero+1;  % it mapped to
            curr_rr = ratio_track_new(curr_i, ind_nonzero(tt));
            %             if curr_j~=ind_nonzero(tt)
            %                 curr_rr = dis_bin_input_new(curr_i)/(curr_j-ind_nonzero(tt));
            %             else
            %                 curr_rr = dis_bin_input_new(curr_i)/0.5;
            %             end
            
            %             tmp = curr_j-ind_nonzero(tt)+1;
            %             if dis_bin_input(curr_i)~=0 && tmp~=0
            %             curr_index = ratio_track(dis_bin_input(curr_i),...
            %                 tmp);
            %             elseif tmp~=0
            %                 curr_index = ratio_track_0(tmp);
            %                 %[qq curr_index] = min(abs(ratio_r-dis_bin_input(curr_i))/(curr_j-ind_nonzero(tt)));
            %             elseif dis_bin_input(curr_i)~=0
            %                 curr_index = ratio_track_1(dis_bin_input(curr_i));
            %             else
            %                 curr_index = 9;
            %                % [qq curr_index] = min(abs(ratio_r-dis_bin_input(curr_i)/0.5));
            %             end
            if pre_rr<0.9220
                pre_qq=1;
            elseif pre_rr>1.085
                pre_qq=3;
            else
                pre_qq=2;
            end
            
            if curr_rr<0.9220
                curr_qq=1;
            elseif curr_rr>1.085
                curr_qq=3;
            else
                curr_qq=2;
            end
            % density_rr = density_mixture_gaussian(log(pre_rr),mix_mu1,mix_std1,mix_mu2,mix_std2,mix_p1,mix_p2);
            %  bbb1= group_transition(curr_qq,pre_qq)*density_rr;
            density_rr =  transition_prob_track(curr_i,curr_j-s+1);
            bbb1 = group_transition(curr_qq,pre_qq)*density_rr;
            
            if bbb1==0
                pp(tt)=0;
            else
                for t=1:length(ind_nonzero)
                    %                     if curr_j~=ind_nonzero(t)
                    %                         curr_rr2=dis_bin_input_new(curr_i)/(curr_j-ind_nonzero(t));
                    %                     else
                    %                         curr_rr2=dis_bin_input_new(curr_i)/0.5;
                    %                     end
                    curr_rr2 = ratio_track_new(curr_i,ind_nonzero(t));
                    
                    %                     tmp = curr_j-ind_nonzero(t)+1;
                    %                     if dis_bin_input(curr_i)~=0 && tmp~=0
                    %                         curr_index2 = ratio_track(dis_bin_input(curr_i),tmp);
                    %                     elseif tmp~=0
                    %                         curr_index2 = ratio_track_0(tmp);
                    %                         %[qq curr_index2] = min(abs(ratio_r-dis_bin_input(curr_i))/(curr_j-ind_nonzero(t)));
                    %                     elseif dis_bin_input(curr_i)~=0
                    %                         curr_index2 = ratio_track_1(dis_bin_input(curr_i));
                    %                     else
                    %                         curr_index2 = 9;
                    %                         %[qq curr_index2] = min(abs(ratio_r-dis_bin_input(curr_i)/0.5));
                    %                     end
                    if curr_rr2<0.9220
                        curr_qq2=1;
                    elseif curr_rr2>1.085
                        curr_qq2=3;
                    else
                        curr_qq2=2;
                    end
                    
                    bbb2= group_transition(curr_qq2,pre_qq)*density_rr;
                    
                    
                    if bbb2~=0
                        ssum2(tt)=ssum2(tt)+exp(log_f{curr_i,curr_j}(ind_nonzero(t))-log_f{curr_i,curr_j}(ind_nonzero(tt)))...
                            *bbb2/bbb1;
                    end;
                end;
                pp(tt)=1/ssum2(tt);
            end;
        end;
        
        candidate=find(pp~=0);
        ssss=0;
        for s=1:size(candidate,1)
            ssss=ssss+pp(candidate(s));
            if u<=ssss
                %  sample_seq(curr_i-1,nn)=
                sample_seq(curr_i-1,nn)=curr_j-ind_nonzero(candidate(s))+1;
                break
                
            end
        end
        s=curr_j-ind_nonzero(candidate(s))+1;
        
        pre_rr = ratio_track_new(curr_i,curr_j-s+1);
        %         if curr_j~=s
        %             pre_rr = dis_bin_input_new(curr_i)/(curr_j-s);
        %         else
        %             pre_rr = dis_bin_input_new(curr_i)/0.5;
        %         end
        %
        %         tmp = curr_j-s+1;
        %         if dis_bin_input_new(curr_i)~=0 && tmp~=0
        %             pre_index = ratio_track(dis_bin_input_new(curr_i),tmp);
        %         elseif tmp~=0
        %             pre_index = ratio_track_0(tmp);
        %             %[qq pre_index] = min(abs(ratio_r-dis_bin_input(curr_i)/(curr_j-s)));
        %         elseif dis_bin_input(curr_i)~=0
        %             pre_index = ratio_track_1(dis_bin_input(curr_i));
        %         else
        %             pre_index = 9;
        %             %[qq pre_index] = min(abs(ratio_r-dis_bin_input(curr_i)/0.5));
        %         end
        curr_i=curr_i-1;
        curr_j=s;
    end
end

for nn=1:size(sample_ali,2)
    for i=1:length(input_scaled_new)
        if sample_seq(i,nn)>0
            sample_ali(i,nn)=discrete_target(sample_seq(i,nn));
        end
    end
end
fprintf('Backward sampling is done.\t\t')
toc