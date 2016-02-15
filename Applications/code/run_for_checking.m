function run_for_checking(coreNAME_woL, dataL, numInput, r_date, exactStart, target_est, varargin)

if nargin < 7
    dataNormalize = 1;
else
    dataNormalize = varargin{1};
end

% if exist('r_date', 'var') == 0
%     r_date = date;
% end
% if exist('exactStart', 'var') == 0
%     exactStart = 4;
% end
% if exist('target_est', 'var') == 0
%     target_est = 'linear';
% end
    
tic
disp('Checking...')

coreNAME = coreNAME_woL;

resultsFolder = ['../../Results/', r_date, '/'];
updateFileN = [coreNAME '_iter'];

preD = dir([resultsFolder updateFileN '*']);
preDN = zeros(length(preD),1);
for i = 1 : length(preD)
    preDN(i) = sscanf(preD(i).name,  [updateFileN '%d']);
end

bb = max(preDN)+1;

resultD = dir([resultsFolder coreNAME '_input*_iter' num2str(bb) '.mat']);
resultDN = zeros(length(resultD),1);
for i = 1 : length(resultD)
   resultDN(i) = sscanf(resultD(i).name, [coreNAME '_input%d_iter' num2str(bb)]);
end

for i = 1 : numInput
    if isempty(find(resultDN == i, 1))
        fid = fopen([coreNAME '_' r_date '_error.txt'], 'a+');
        fprintf(fid, 'No results of input %d\n', i);
        disp(['No results of input ' num2str(i)])
        fclose(fid);
    end
end


