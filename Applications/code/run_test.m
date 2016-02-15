% for totalI = 1 : 3
%     run_for_communication('LR04_5', 50, '18-Apr-2015', 0)
%     for inputNN = 1 : 5
%         run_for_each_input('LR04_5', 50, inputNN, '18-Apr-2015', 0)
%     end
% end
clear all
run_for_communication('S300_7', 100, '22-Apr-2015', 0, 'step')
for totalI = 1 : 20
    for inputNN = 1 : 20
        run_for_each_input('S300_7', 100, inputNN, '22-Apr-2015', 0, 'step')
    end
    run_for_communication('S300_7', 100, '22-Apr-2015', 0, 'step')
end

clear all
run_for_communication('S300_7', 100, '22-Apr-2015', 0, 'linear')
for totalI = 1 : 20
    for inputNN = 1 : 20
        run_for_each_input('S300_7', 100, inputNN, '22-Apr-2015', 0, 'linear')
    end
    run_for_communication('S300_7', 100, '22-Apr-2015', 0, 'linear')
end

%%
% inputList = dir('LR04_2_250_IC0_input*_bb4.mat')
% 
% figure
% hold on
% for i = 1 : length(inputList)
%     load(inputList(i).name)
% %     totalsample{i} = sample_seq;
%     plot(sample_ali(:,1), input_scaled_new(:,2))
% end
    