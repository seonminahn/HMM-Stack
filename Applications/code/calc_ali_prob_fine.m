%returns several statistics for each alignment provided
%note: input ali_seq requires to enter alignment with aligned bin index
%infomation.
% - log_prob: log likeligood of the alignment
% - ali_ratio: mapping rate value at each input data point
% - ratio_index_seq: mapping rate index info at each input data point
% - log_residual_seq: log of pdf of residual at each input data point
% - log_rate_change_seq: log of transition probability at each input data
% point

function [log_prob ali_ratio log_residual_seq log_rate_change_seq] = calc_ali_prob_fine(ali_seq,input_scaled_new,pi_0,pi_nomatch,dis_bin_input_new,ratio_track,log_grid,ratio_track_new,transition_prob_track,group_transition)


ali_ratio = zeros(length(input_scaled_new)-1,1);
log_residual_seq = zeros(length(input_scaled_new),1);
log_rate_change_seq = zeros(length(input_scaled_new)-1,1);

have_align = find(ali_seq~=0);
align_start = have_align(1);
align_end = have_align(end);
log_prob = 0;
if (align_start-1)+(size(log_grid,2)-ali_seq(align_end))~=0
    log_prob = log(pi_nomatch)*((align_start-1)+(size(log_grid,2)-ali_seq(align_end)))+...
        log_grid(align_start,ali_seq(align_start));
end


log_residual_seq(align_start) = log_grid(align_start,ali_seq(align_start));
i = align_start+1;
tmp = ali_seq(i)-ali_seq(i-1)+1;

pre_index = ratio_track(dis_bin_input_new(i),tmp);
% if tmp~=0 && dis_bin_input_new(i)~=0
%     pre_index = ratio_track(dis_bin_input_new(i),tmp);
%     %elseif tmp~=0
%     %   pre_index = ratio_track_0(tmp);
% elseif dis_bin_input_new(i)~=0
%     pre_index = ratio_track_1(dis_bin_input_new(i));
% else
%     pre_index = 9;
% end
pre_rr = ratio_track_new(i,ali_seq(i)-ali_seq(i-1)+1);
% if ali_seq(i)-ali_seq(i-1)~=0
%     pre_rr = dis_bin_input_new(i)/(ali_seq(i)-ali_seq(i-1));
% else
%     pre_rr = dis_bin_input_new(i)/0.5;
% end
%pre_index = ratio_track(dis_bin_input_new(i),ali_seq(i)-ali_seq(i-1)+1);
ali_ratio(i-1) = pre_rr;
log_rate_change_seq(i-1) = log(pi_0(pre_index));
log_prob = log_prob + log(pi_0(pre_index))+log_grid(i,ali_seq(i));

log_residual_seq(i)=log_grid(i,ali_seq(i));

for i=align_start+2:align_end
%     tmp = ali_seq(i)-ali_seq(i-1);
%     if tmp~=0 && dis_bin_input_new(i)~=0
%         curr_index = ratio_track(dis_bin_input_new(i),tmp);
%         %elseif tmp~=0
%         %   curr_index = ratio_track_0(tmp);
%     elseif dis_bin_input_new(i)~=0
%         curr_index = ratio_track_1(dis_bin_input_new(i));
%     else
%         curr_index = 9;
%     end
    curr_rr = ratio_track_new(i,ali_seq(i)-ali_seq(i-1)+1);
%     if ali_seq(i)~=ali_seq(i-1)
%         curr_rr = dis_bin_input_new(i)/(ali_seq(i)-ali_seq(i-1));
%     else
%         curr_rr = dis_bin_input_new(i)/0.5;
%     end
    
    if curr_rr<0.9220
        curr_qq=1;
    elseif curr_rr>1.085
        curr_qq=3;
    else
        curr_qq=2;
    end
    
    
    if pre_rr<0.9220
        pre_qq=1;
    elseif pre_rr>1.085
        pre_qq=3;
    else
        pre_qq=2;
    end
    density_rr = transition_prob_track(i,ali_seq(i)-ali_seq(i-1)+1);
%     density_rr = density_mixture_gaussian(log(curr_rr),mix_mu1,mix_std1,mix_mu2,mix_std2,mix_p1,mix_p2);
     log_A_new= log(group_transition(pre_qq,curr_qq)*density_rr);
    
    
    %curr_index = ratio_track(dis_bin_input_new(i),ali_seq(i)-ali_seq(i-1)+1);
    log_prob = log_prob + log_A_new...
        +log_grid(i,ali_seq(i));
    log_rate_change_seq(i-1) = log_A_new;
    pre_rr = curr_rr;
    ali_ratio(i-1) = pre_rr;
    log_residual_seq(i)=log_grid(i,ali_seq(i));
    
end
% 
% ali_ratio = zeros(length(ratio_index_seq),1);
% for i=1:length(ratio_index_seq)
%     if ratio_index_seq(i)~=0
%         ali_ratio(i)=ratio_r(ratio_index_seq(i));
%     end
% end