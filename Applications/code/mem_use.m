function mem = mem_use()
% get approximate memory use in GB on Linux

pid=feature('getpid');  % process ID
[~, output]=system(['cat /proc/', num2str(pid), '/statm']);
m  = regexp(output, ('\d+'), 'match');

mem = str2num(m{1}) * 4 /1024/1024;