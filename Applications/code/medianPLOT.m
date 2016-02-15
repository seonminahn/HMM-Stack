currentF = '/Users/seonminahn/Dropbox/Brown_University/ALIGN/ProfileHMM/code_script_run';
resultF = '/Users/seonminahn/Documents/Brown_University/2015/ALIGN/ProfileHMM_results/26-May-2015/bb7';

for medianFI = [3 5 6 7 16 17 76 81 82 87]
    
    cd(resultF)
    
    fileNmedian = ['LR04plusNew_6000_IC1_step_norm0_input',  num2str(medianFI), '_bb7_median.mat'];
    if exist(fileNmedian, 'file') == 2
        medianFI
        
        load(['LR04plusNew_6000_IC1_step_norm0_input',  num2str(medianFI), '_bb7.mat'])
        load(fileNmedian)
        
        cd(currentF)
        plot_align
        
        cd(resultF)
        savefig(['LR04plusNew_6000_IC1_step_norm0_input',  num2str(medianFI), '_bb7.fig'])
        export_fig(['LR04plusNew_6000_IC1_step_norm0_input',  num2str(medianFI), '_bb7.pdf'])
        
        
        h2 = figure;
        hold on
        plot(median_ali, input_scaled_new(:,2), 'bo-')
        stairs(stairsX, total_mu(unique_bin_mean_index),'r-')
        % stairs(total_target(unique_bin_mean_index), total_mu(unique_bin_mean_index),'ro')
        set(gca,'Ydir','reverse')
        legend(['Core: ' inputNAME_name{inputN}], 'New Stack')
        title('Sequences after alignment');
        xlabel('Stack Age (kyr)')
        ylabel('\delta^{18}O')
        axis tight
        savefig(['LR04plusNew_6000_IC1_step_norm0_input',  num2str(medianFI), '_bb7_median.fig'])
        export_fig(['LR04plusNew_6000_IC1_step_norm0_input',  num2str(medianFI), '_bb7_median.pdf'])
        
%         pause
        close(h1)
        close(h2)
    end
    
    
end


%%
% Time Collect
resultF = '/Users/seonminahn/Documents/Brown_University/2015/ALIGN/ProfileHMM_results/26-May-2015/bb7/data';

for medianFI = 1:196
    
    cd(resultF)
    
    fileNmedian = ['LR04plusNew_6000_IC1_step_norm0_input',  num2str(medianFI), '_bb7_median.mat'];
    if exist(fileNmedian, 'file') == 2
        load(fileNmedian)
        medianT(medianFI,1) = calT;
        medianFI
    end
end

medianT(:,2) = medianT(:,1)/60;
medianT(:,3) = medianT(:,2)/60;


%% Find name
for i = 1 : 196
   if strcmp(inputNAME_name{i}, 'MD02-2575\_3') 
        i
   end
end

%% input length
for i = 1 : 196
    load(['LR04plusNew_6000_IC1_step_norm0_input', num2str(i), '_bb7_updateD.mat'])
    l(i,1) = input_scaled_new(end,1) - input_scaled_new(1,1);
    inputS(i,1) = input_scaled_new(1,1);
    inputE(i,1) = input_scaled_new(end,1);
end

%% 
(find(l>1000));
for i = li'
    figFile = ['LR04plusNew_6000_IC1_step_norm0_input', num2str(i), '_bb7_median.fig'];
    if exist(figFile, 'file')
    open(figFile)
    pause
%         close(figure(2))
    end
end