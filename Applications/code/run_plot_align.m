clear 
fList = dir;

for i = 4 : length(fList)
    fildN(i) = sscanf(fList(i).name, 'LR04_5325_IC1_step_norm0_input%d');
end
fildN(1:3) = [];

resultF = '/Users/seonminahn/Documents/Results/30-Sep-2015-median/';

%%
for fI = 4 : length(fList)
   clearvars -except resultF fList fI
   load([resultF fList(fI).name])
   plot_align
   
   figN = strrep(fList(fI).name, 'medain.mat', 'median');
   
   figure(1); savefig(h1, [resultF figN '_1.fig']); export_fig([resultF figN '_1.pdf']);
   figure(2); savefig(h2, [resultF figN '_2.fig']); export_fig([resultF figN '_2.pdf']);
   close(h1)
   close(h2)
end