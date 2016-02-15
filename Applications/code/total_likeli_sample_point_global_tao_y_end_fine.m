%calculate the total likelihood of the samples to maximize the expectation
%in EM

function total_likeli_prob = total_likeli_sample_point_global_tao_y_end_fine(sample_seq,pi_0,log_pi_begin,log_tao_x_begin,log_tao_y_begin,log_tao_x_end,log_delta_match,log_delta_end,log_grid,ratio_track,dis_bin_input_new,log_input_gap_raw,aaa,transition_prob_track,ratio_track_new,group_transition)
log_tao_y_end = [log(aaa) log(1-aaa)];

% Evaluate the likelihood with a new parameter
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
    
    log_prob = log_grid(align_start,ali_seq(align_start)) + log(pi_0(pre_index))+log_grid(i,ali_seq(i))+log_delta_match(1);
    
    for i=align_start+2:align_end
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
        log_prob = log_prob + log_delta_match(2)+log_delta_end(1);
    elseif align_end==size(sample_seq,1);
        log_prob = log_prob + log_delta_match(2)+log_delta_end(3)+log_tao_y_end(1)*(size(log_grid,2)-ali_seq(end)-1)+log_tao_y_end(2);
    else
        log_prob = log_prob + log_delta_match(2)+log_delta_end(2)+log_tao_x_end(1)*(size(log_grid,1)-align_end-1)+log_tao_x_end(2);
        for j=align_end+1:length(log_input_gap_raw)
            log_prob = log_prob +log_input_gap_raw(j);
        end
    end
    likeli_sample(s)=log_prob;    
end

total_likeli_prob = -sum(likeli_sample);