%calculate the total likelihood of the samples to maximize the expectation
%in EM

function total_likeli_prob = total_likeli_sample_point_global_min_delta_r_fine(sample_seq,pi_0,log_pi_begin,log_tao_x_begin,log_tao_y_begin,log_tao_x_end,log_tao_y_end,log_grid,ratio_track,dis_bin_input_new,log_input_gap_raw,aaa,transition_prob_track,ratio_track_new,group_transition)
delta_r =[aaa(1) aaa(2) 1-aaa(1)-aaa(2)];
log_delta_r = log(delta_r);
%log_A = log(A_enlarge);
likeli_sample = zeros(size(sample_seq,2),1);

for s=1:size(sample_seq,2)
    
    ali_seq = sample_seq(:,s);
    have_align = find(ali_seq~=0);
    align_start = have_align(1);
    align_end = have_align(end);
    i = align_start+1;
    tmp = ali_seq(i)-ali_seq(i-1)+1;
    
    pre_index = ratio_track(dis_bin_input_new(i),tmp);
    pre_rr = ratio_track_new(i,tmp);
%     if dis_bin_input_new(i)~=0 && tmp~=0
%         pre_index = ratio_track(dis_bin_input_new(i),tmp);
%         %elseif tmp~=0
%         %   pre_index = ratio_track_0(tmp);
%     elseif dis_bin_input_new(i)~=0
%         pre_index = ratio_track_1(dis_bin_input_new(i));
%     else
%         pre_index=9;
%         
%     end
%     
%     if ali_seq(i)-ali_seq(i-1)~=0
%         pre_rr = dis_bin_input_new(i)/(ali_seq(i)-ali_seq(i-1));
%     else
%         pre_rr = dis_bin_input_new(i)/0.5;
%     end
    %pre_index = ratio_track(dis_bin_input_new(i),ali_seq(i)-ali_seq(i-1)+1);
    log_prob = log_grid(align_start,ali_seq(align_start)) + log(pi_0(pre_index))+log_grid(i,ali_seq(i));
    for i=align_start+2:align_end
        %         tmp = ali_seq(i)-ali_seq(i-1)+1;
        %         if dis_bin_input_new(i)~=0 && tmp~=0
        %             curr_index = ratio_track(dis_bin_input_new(i),tmp);
        %         elseif tmp~=0
        %             curr_index=ratio_track_0(tmp);
        %         elseif dis_bin_input_new(i)~=0
        %             curr_index=ratio_track_1(dis_bin_input_new(i));
        %         else
        %             curr_index=9;
        %
        %         end
%         if ali_seq(i)-ali_seq(i-1)~=0
%             curr_rr = dis_bin_input_new(i)/(ali_seq(i)-ali_seq(i-1));
%         else
%             curr_rr = dis_bin_input_new(i)/0.5;
%         end
%         
        curr_rr = ratio_track_new(i,ali_seq(i)-ali_seq(i-1)+1);
        
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
        density_rr = transition_prob_track(i,ali_seq(i)-ali_seq(i-1)+1);
%         density_rr = density_mixture_gaussian(log(curr_rr),mix_mu1,mix_std1,mix_mu2,mix_std2,mix_p1,mix_p2);
%         %ssum= group_transition(qq1,qq2)*density_rr;
        
        %curr_index = ratio_track(dis_bin_input_new(i),ali_seq(i)-ali_seq(i-1)+1);
        log_prob = log_prob + log(group_transition(pre_qq,curr_qq)*density_rr)+log_grid(i,ali_seq(i))+log_delta_match(1);
        pre_rr = curr_rr;
        
    end
    
    
    if align_start==1 && ali_seq(1)==1
        log_prob = log_prob+log_pi_begin(2);
    elseif align_start~=1
        log_prob = log_prob + log_pi_begin(1)+log_tao_x_begin(1)*(align_start-2)+log_tao_x_begin(2);
        for j=1:align_start-1
            log_prob = log_prob +log_input_gap_raw(j);
        end
    else
        log_prob = log_prob + log_pi_begin(3)+log_tao_y_begin(1)*(ali_seq(1)-2)+log_tao_y_begin(2);
    end
    
    if align_end==size(sample_seq,1) && ali_seq(size(sample_seq,1))==size(log_grid,2)
        log_prob = log_prob + log_delta_r(3);
    elseif align_end==size(sample_seq,1);
        log_prob = log_prob + log_delta_r(3)+log_tao_y_end(1)*(size(log_grid,2)-ali_seq(end)...
            -1)+log_tao_y_end(2);
    else
        log_prob = log_prob + log_delta_r(2)+log_tao_x_end(1)*(size(log_grid,1)-align_end-1)...
            +log_tao_x_end(2);
        for j=align_end+1:length(log_input_gap_raw)
            log_prob = log_prob +log_input_gap_raw(j);
        end
    end
    likeli_sample(s)=log_prob;
    
end

total_likeli_prob = -sum(likeli_sample);
