% %%
% clear
% selectInput = 32;
% 
% currentF = pwd;
% resultF = '/Volumes/ANNEVIEW/ProfileHMM_results/26-May-2015/bb7';
% 
% cd(resultF)
% load(['LR04plusNew_6000_IC1_step_norm0_input', num2str(selectInput), '_bb7.mat'])
% % load(['LR04plusNew_6000_IC1_step_norm0_input', num2str(selectInput), '_bb7_median.mat'])
% cd(currentF)
% %%
% % cd /Users/seonminahn/Dropbox/Brown_University/ALIGN/HMM_Match/HMM_Match_201411_original
% % cd /Users/seonminahn/Dropbox/Brown_University/ALIGN/ProfileHMM/code_script_run
% %obtain median estimator
% tic
% calc_median_ali_dp_fine_half;
% toc
% cd(resultF)
% save(['LR04plusNew_6000_IC1_step_norm0_input', num2str(selectInput), '_bb7_median.mat'], 'median_seq', 'median_ali')
% cd(currentF)

%%
%calculate upper and lower limit
disp('Drawing age estimates...')

[upper_95 lower_95] = upper_lower_limit(input_scaled_new, sample_ali);


if exist('meanShift') == 0
    meanShift = 0;
end

for i = 1 : length(inputNAME)
   str = inputNAME{i};
   str = strrep(str, '../NewBenthicdata/NewRecords/', '');
   str = strrep(str, '../NewBenthicdata/LR04/', '');
   str = strrep(str, '_LR04age', '');
   str = strrep(str, '.txt', '');
   str = strrep(str, '_', '\_');
   str = strrep(str, '../data/', '');
   inputNAME_name{i,1} = str;
end



h1 = figure;

unique_bin_mean_index = unique(bin_mean_index);
stairsX = [total_target(unique_bin_mean_index(1),1); total_target(unique_bin_mean_index(1:end-1),1)+diff(total_target(unique_bin_mean_index))/2];

% subplot(2,2,1)
% hold on
% plot(input_scaled_new(:,1),input_scaled_new(:,2),'bo-')
% stairs(stairsX, total_mu(unique_bin_mean_index),'r-')
% % stairs(total_target(unique_bin_mean_index), total_mu(unique_bin_mean_index),'ro')
% set(gca,'Ydir','reverse')
% l = legend(['Core: ' inputNAME_name{inputN}], 'HMM\_LR04');
% title('Sequences before alignment');
% % xlabel('Stack Age (kyr)')
% set(gca, 'XTick', [])
% ylabel(['\delta^{18}O (' char(8240) ')'])
% axis tight

% subplot(2,2,3)
subplot(3,1,1)
hold on
plot(median_ali, input_scaled_new(:,2)+meanShift, 'bo-')
stairs(stairsX, total_mu(unique_bin_mean_index),'r-')
% stairs(total_target(unique_bin_mean_index), total_mu(unique_bin_mean_index),'ro')
set(gca,'Ydir','reverse')
% legend(['Core: ' inputNAME_name{inputN}], 'New Stack')
% title('Sequences after alignment');
% xlabel('Age (ka)')
ylabel(['\delta^{18}O (' char(8240) ')'])
axis tight
% l = legend(['Core: ' inputNAME_name{inputN}], 'HMM\_LR04');
l = legend(['Core: ' inputNAME_name{inputN}], 'HMM stack');
l.Box = 'off';
g = gca;
g.Box = 'on';


% subplot(2,2,2)
subplot(3,1,2)
ali = median_ali;
aa = find(ali~=0);
align_start=aa(1);
align_end=aa(end);
ali_new=ali;
for i=align_end+1:length(ali)
    ali_new(i)=ali(align_end);
end
plot(ali_new,ali-ali,'r',ali_new,upper_95-ali_new,'b',ali_new,lower_95-ali_new,'g')
ll = legend('median','upper limit', 'lower limit');
ll.Box = 'off';
%title('point-version -- 95% confidence limit of aligned age (normalization) at input')
% xlabel('Age (ka)')
% title('Confidence Band in Age (kyr)')
title('95% interval of age estimates')
ylabel('Age (kyr)')
axis tight

% subplot(2,2,4)
subplot(3,1,3)
ali_seq = median_seq;
[log_prob ali_ratio log_residual_seq log_rate_change_seq] = ...
    calc_ali_prob_fine(ali_seq,input_scaled_new,pi_0,pi_nomatch,dis_bin_input_new,...
    ratio_track,log_grid,ratio_track_new,transition_prob_track,group_transition);
ali = zeros(length(input_scaled_new),1);
for i=1:length(input_scaled_new)
    if ali_seq(i)~=0
        ali(i)=discrete_target(ali_seq(i));
    end
end
vv=find(ali(2:end)~=0);
semilogy(ali(vv),ali_ratio(vv),'bo-')
xlabel('Age (ka)')
ylabel('Rate')
title('Relative Accumulation Rates')
axis tight
set(gca, 'YTickLabel',{'0','0.5','1','2','3','4'}, 'YTick',[0 0.5 1 2 3 4])
currentA = axis;
currentA(3:4) = [0 4];
axis(currentA)


% h2 = figure;
% hold on
% plot(median_ali, input_scaled_new(:,2)+meanShift, 'bo-')
% stairs(stairsX, total_mu(unique_bin_mean_index),'r-')
% % stairs(total_target(unique_bin_mean_index), total_mu(unique_bin_mean_index),'ro')
% set(gca,'Ydir','reverse')
% % legend(['Core: ' inputNAME_name{inputN}], 'New Stack')
% title('Sequences after alignment');
% xlabel('Stack Age (kyr)')
% ylabel('\delta^{18}O')
% axis tight

% a1 = annotation('textbox', [0.1 0.91 0.03 0.05], 'String', 'a', 'Fontsize', 24, 'FontWeight', 'bold', 'LineStyle', 'none');
% a2 = annotation('textbox', [0.1 0.435 0.03 0.05], 'String', 'c', 'Fontsize', 24, 'FontWeight', 'bold', 'LineStyle', 'none');
% a3 = annotation('textbox', [0.5424 0.91 0.03 0.05], 'String', 'b', 'Fontsize', 24, 'FontWeight', 'bold', 'LineStyle', 'none');
% a4 = annotation('textbox', [0.5424 0.435 0.03 0.05], 'String', 'd', 'Fontsize', 24, 'FontWeight', 'bold', 'LineStyle', 'none');
% ll = l.Position;
% l.Position = [0.37 0.5 ll(3) ll(4)];

a1 = annotation('textbox', [0.1 0.9 0.03 0.05], 'String', 'a', 'Fontsize', 24, 'FontWeight', 'bold', 'LineStyle', 'none');
a2 = annotation('textbox', [0.1 0.3 0.03 0.05], 'String', 'c', 'Fontsize', 24, 'FontWeight', 'bold', 'LineStyle', 'none');
a3 = annotation('textbox', [0.1 0.6 0.03 0.05], 'String', 'b', 'Fontsize', 24, 'FontWeight', 'bold', 'LineStyle', 'none');

clear h1
