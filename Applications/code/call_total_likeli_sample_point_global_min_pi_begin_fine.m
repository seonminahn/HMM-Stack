function xx = call_total_likeli_sample_point_global_min_pi_begin_fine(sample_seq,pi_0,log_tao_x_begin,log_tao_y_begin,log_tao_x_end,log_tao_y_end,log_delta_match,log_delta_end,ratio_track,dis_bin_input_new,log_grid,log_input_gap_raw ,transition_prob_track,ratio_track_new,group_transition)

fhandle = @(aaa) total_likeli_sample_point_global_min_pi_begin_fine(sample_seq,pi_0,log_tao_x_begin,log_tao_y_begin,log_tao_x_end,log_tao_y_end,log_delta_match,log_delta_end,ratio_track,dis_bin_input_new,log_grid,log_input_gap_raw,aaa,transition_prob_track,ratio_track_new,group_transition);
options = optimoptions('fmincon', 'display', 'notify-detailed');
%SA% xx = fmincon(fhandle,[0.45; 0.45],[1 1],0.999,[],[],[0.001;0.001],[0.999;0.999],[], options);
xx = fmincon(fhandle,[0.33; 0.33],[1 1],0.999,[],[],[0.001;0.001],[0.999;0.999],[],options);

end
