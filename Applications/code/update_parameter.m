%%%SA%%%
fhandle = @(new_meanShift) total_likeli_sample_point_global_min_fine(sample_seq,pi_0,log_pi_begin,log_tao_x_begin,log_tao_y_begin,log_tao_x_end,log_tao_y_end,log_delta_match,log_delta_end,ratio_track,dis_bin_input_new,mu,sigma,delta_grid,input_scaled_new,mean_tar,transition_prob_track,ratio_track_new,group_transition,new_meanShift);
[new_meanShift, ~] = fminbnd(fhandle,-2,2);
mu = mu - new_meanShift;
mean_mu = mean_mu - new_meanShift;
meanShift = meanShift + new_meanShift;

log_grid = log(normpdf(delta_grid, ones(length(input_scaled_new),1)*mu', ones(length(input_scaled_new),1)*sigma'));

log_input_gap_raw = zeros(length(input_scaled_new),1);
for i=1:length(input_scaled_new)
    tempNormpdf = normpdf(input_scaled_new(i,2)-mean_tar, mean_mu, mean_sigma);
    log_input_gap_raw(i) = log(tempNormpdf);
end


xx = call_total_likeli_sample_point_global_min_pi_begin_fine(sample_seq,pi_0,log_tao_x_begin,log_tao_y_begin,log_tao_x_end,log_tao_y_end,log_delta_match,log_delta_end,ratio_track,dis_bin_input_new,log_grid,log_input_gap_raw,transition_prob_track,ratio_track_new,group_transition);
pi_begin = [xx(1) xx(2) 1-xx(1)-xx(2)];
log_pi_begin = log(pi_begin);

bbb = call_total_likeli_sample_point_global_tao_x_begin_fine(sample_seq,pi_0,log_pi_begin,log_tao_y_begin,log_tao_x_end,log_tao_y_end,log_delta_match,log_delta_end,log_grid,ratio_track,dis_bin_input_new,log_input_gap_raw,transition_prob_track,ratio_track_new,group_transition);
tao_x_begin =[bbb 1-bbb];
log_tao_x_begin = log(tao_x_begin);

bbb = call_total_likeli_sample_point_global_tao_y_begin_fine(sample_seq,pi_0,log_pi_begin,log_tao_x_begin,log_tao_x_end,log_tao_y_end,log_delta_match,log_delta_end,log_grid,ratio_track,dis_bin_input_new,log_input_gap_raw,transition_prob_track,ratio_track_new,group_transition);
tao_y_begin =[bbb 1-bbb];
log_tao_y_begin = log(tao_y_begin);

bbb = call_total_likeli_sample_point_global_tao_x_end_fine(sample_seq,pi_0,log_pi_begin,log_tao_x_begin,log_tao_y_begin,log_tao_y_end,log_delta_match,log_delta_end,log_grid,ratio_track,dis_bin_input_new,log_input_gap_raw,transition_prob_track,ratio_track_new,group_transition);
tao_x_end =[bbb 1-bbb];
log_tao_x_end = log(tao_x_end);

bbb = call_total_likeli_sample_point_global_tao_y_end_fine(sample_seq,pi_0,log_pi_begin,log_tao_x_begin,log_tao_y_begin,log_tao_x_end,log_delta_match,log_delta_end,log_grid,ratio_track,dis_bin_input_new,log_input_gap_raw,transition_prob_track,ratio_track_new,group_transition);
tao_y_end =[bbb 1-bbb];
log_tao_y_end = log(tao_y_end);

% xx = call_total_likeli_sample_point_global_min_delta_r_fine(sample_seq,pi_0,log_pi_begin,log_tao_x_begin,log_tao_y_begin,log_tao_x_end,log_tao_y_end,log_grid,ratio_track,dis_bin_input_new,log_input_gap_raw,transition_prob_track,ratio_track_new,group_transition);
% delta_r = [xx(1) xx(2) 1-xx(1)-xx(2)];
% log_delta_r = log(delta_r);

fhandle = @(aaa) total_likeli_sample_point_global_min_delta_end(sample_seq,pi_0,log_pi_begin,log_tao_x_begin,log_tao_y_begin,log_tao_x_end,log_tao_y_end,log_delta_match,log_grid,ratio_track,dis_bin_input_new,log_input_gap_raw,aaa,transition_prob_track,ratio_track_new,group_transition);
options = optimoptions('fmincon', 'display', 'notify-detailed');
xx = fmincon(fhandle,[0.33; 0.33],[1 1],0.999,[],[],[0.001;0.001],[0.999;0.999],[],options);
delta_end = [xx(1) xx(2) 1-xx(1)-xx(2)];
log_delta_end = log(delta_end);

fhandle = @(aaa) total_likeli_sample_point_global_min_delta_match(sample_seq,pi_0,log_pi_begin,log_tao_x_begin,log_tao_y_begin,log_tao_x_end,log_tao_y_end,log_delta_end,log_grid,ratio_track,dis_bin_input_new,log_input_gap_raw,aaa,transition_prob_track,ratio_track_new,group_transition);
[bbb,fval] = fminbnd(fhandle,0.001 ,0.999);
delta_match =[bbb 1-bbb];
log_delta_match = log(delta_match);

%         fhandle = @(aaa) total_likeli_sample_point_global_min_pi_begin_fine(sample_seq,pi_0,log_tao_x_begin,log_tao_y_begin,log_tao_x_end,log_tao_y_end,log_delta_r,ratio_track,dis_bin_input_new,log_grid,log_input_gap_raw,aaa,transition_prob_track,ratio_track_new,group_transition);
%         xx = fmincon(fhandle,[0.45; 0.45],[1 1],0.999,[],[],[0.001;0.001],[0.999;0.999]);
%         pi_begin = [xx(1) xx(2) 1-xx(1)-xx(2)];
%         log_pi_begin = log(pi_begin);
%
%         fhandle = @(aaa) total_likeli_sample_point_global_tao_x_begin_fine(sample_seq,pi_0,log_pi_begin,log_tao_y_begin,log_tao_x_end,log_tao_y_end,log_delta_r,log_grid,ratio_track,dis_bin_input_new,log_input_gap_raw,aaa,transition_prob_track,ratio_track_new,group_transition);
%         [bbb,fval] = fminbnd(fhandle,0.001 ,0.999);
%         tao_x_begin =[bbb 1-bbb];
%         log_tao_x_begin = log(tao_x_begin);
%
%         fhandle = @ (aaa) total_likeli_sample_point_global_tao_y_begin_fine(sample_seq,pi_0,log_pi_begin,log_tao_x_begin,log_tao_x_end,log_tao_y_end,log_delta_r,log_grid,ratio_track,dis_bin_input_new,log_input_gap_raw,aaa,transition_prob_track,ratio_track_new,group_transition);
%         [bbb, fval] = fminbnd(fhandle,0.001 ,0.999);
%         tao_y_begin =[bbb 1-bbb];
%         log_tao_y_begin = log(tao_y_begin);
%
%         fhandle = @ (aaa) total_likeli_sample_point_global_tao_x_end_fine(sample_seq,pi_0,log_pi_begin,log_tao_x_begin,log_tao_y_begin,log_tao_y_end,log_delta_r,log_grid,ratio_track,dis_bin_input_new,log_input_gap_raw,aaa,transition_prob_track,ratio_track_new,group_transition);
%         [bbb, fval] = fminbnd(fhandle,0.001 ,0.999);
%         tao_x_end =[bbb 1-bbb];
%         log_tao_x_end = log(tao_x_end);
%
%         fhandle = @ (aaa) total_likeli_sample_point_global_tao_y_end_fine(sample_seq,pi_0,log_pi_begin,log_tao_x_begin,log_tao_y_begin,log_tao_x_end,log_delta_r,log_grid,ratio_track,dis_bin_input_new,log_input_gap_raw,aaa,transition_prob_track,ratio_track_new,group_transition);
%         [bbb, fval] = fminbnd(fhandle,0.001 ,0.999);
%         tao_y_end =[bbb 1-bbb];
%         log_tao_y_end = log(tao_y_end);
%
%         fhandle = @(aaa) total_likeli_sample_point_global_min_delta_r_fine(sample_seq,pi_0,log_pi_begin,log_tao_x_begin,log_tao_y_begin,log_tao_x_end,log_tao_y_end,log_grid,ratio_track,dis_bin_input_new,log_input_gap_raw,aaa,transition_prob_track,ratio_track_new,group_transition);
%         xx = fmincon(fhandle,[0.45; 0.45],[1 1],0.999,[],[],[0.001;0.001],[0.999;0.999]);
%         delta_r = [xx(1) xx(2) 1-xx(1)-xx(2)];
%         log_delta_r = log(delta_r);
%%%%%%%%

fre_start=zeros(1,rmax);
fre_A = zeros(rmax,rmax);
ratio_seq = cell(size(sample_seq,2),1);
for nn=1:size(sample_seq,2)
    sam = sample_seq(:,nn);
    have_align = find(sam~=0);
    align_start = have_align(1);
    align_end = have_align(end);
    ratio_seq{nn} = zeros(align_end-align_start,1);
    for i=align_start+1:align_end
        curr_rr = ratio_track_new(i,sam(i)-sam(i-1)+1);
        %        if sam(i)~=sam(i-1)
        %            curr_rr = dis_bin_input_new(i)/(sam(i)-sam(i-1));
        %        else
        %            curr_rr = dis_bin_input_new(i)/0.5;
        %        end
        %
        
        %        tmp = sam(i)-sam(i-1);
        %        if dis_bin_input_new(i)~=0 && tmp~=0
        %       nearest_index = ratio_track(dis_bin_input_new(i),tmp);
        %       % elseif tmp~=0
        %        %    nearest_index = ratio_track_0(tmp);
        %            %[qq nearest_index] = min(abs(ratio_r-dis_bin_input(i)/(sam(i)-sam(i-1))));
        %        elseif dis_bin_input_new(i)~=0
        %            nearest_index = ratio_track_1(dis_bin_input_new(i));
        %        else
        %            nearest_index = 9;
        %           % [qq nearest_index] = min(abs(ratio_r-dis_bin_input(i)/0.5));
        %        end
        %nearest_index = ratio_track(dis_bin_input(i),sam(i)-sam(i-1)+1);
        ratio_seq{nn}(i-align_start) = curr_rr;
    end
end
total=size(sample_seq,2);
for s=1:rmax
    for ss=1:size(sample_seq,2)
        temp_seq=ratio_seq{ss};
        if temp_seq(1)==s
            fre_start(s)=fre_start(s)+1;
        end;
        %%%SA%%%
        temp_seq = [];
        % clear temp_seq;
        %%%%%%%%
    end;
    if fre_start(s)==0
        fre_start(s)=1;
    end;
    total=total+1;
end;

% for ss=1:size(sample_seq,2)
%     temp_seq=ratio_seq{ss};
%     for nn=1:length(temp_seq)-1
%         fre_A(temp_seq(nn),temp_seq(nn+1))=fre_A(temp_seq(nn),temp_seq(nn+1))+1;
%     end;
%     clear temp_seq;
% end;

%add pseudo count
% for i=1:rmax
%     for j=1:rmax
%         if fre_A(i,j)==0
%             fre_A(i,j)=1;
%         end;
%     end;
% end;


pi_0=zeros(rmax,1);
for s=1:rmax
    pi_0(s)=fre_start(s)/sum(fre_start);
end;


%%%SA%%%
%update the log_grid after updating sigma
% log_grid=log(normpdf(delta_grid,mu,sigma));
% log_input_gap_raw = zeros(length(input_scaled_new),1);
% for i=1:length(input_scaled_new)
%     log_input_gap_raw(i) = log(normpdf(input_scaled_new(i,2)-mean_tar,mu,sigma));
% end
%%%%%%%%


fprintf('Updating parameters is done.\t\t')
toc


% %% Display new estimates
% fprintf('log likelihood in iteration %d is %f\n',bb, log_fE);
% %         fprintf('The sigma returned in iteration %d is %f\n', bb, sigma);
% %         fprintf('The mu returned in iteration %d is %f\n', bb, mu);
% fprintf('The pi_begin returned in iteration %d is %f %f\n', bb, pi_begin(1), pi_begin(2));
% fprintf('The tao_x_begin returned in iteration %d is %f\n', bb, tao_x_begin(1));
% fprintf('The tao_x_end returned in iteration %d is %f\n', bb, tao_x_end(1));
% fprintf('The tao_y_begin returned in iteration %d is %f\n', bb, tao_y_begin(1));
% fprintf('The tao_y_end returned in iteration %d is %f\n', bb, tao_y_end(1));
% fprintf('The delta_r returned in iteration %d is %f %f\n', bb, delta_r(1), delta_r(2));
