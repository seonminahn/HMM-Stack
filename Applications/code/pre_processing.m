mean_tar = 0;
qqqq=1;

%%%SA%%%
% Sort seq1 in depth 
[~, si] = sort(seq1(:,1));
seq1 = seq1(si,:);

% Choose the first one of duplicated values
% [input(:,1),label]=unique(seq1(:,1));
% input(:,2)=seq1(label,2);
tmp=unique(seq1(:,1));
input=zeros(length(tmp),2);
for i=1:length(tmp)
    % Choose one of duplicated values randomly
    d_input = find(seq1(:,1)==tmp(i));
    chooseI = min(d_input) + ceil(rand*(length(d_input))) - 1;
    input(i,:)= [tmp(i) seq1(chooseI,2)];
    
    % Choose the average of duplicated values
    %     input(i,:)= [tmp(i) mean(seq1(seq1(:,1)==tmp(i),2))];
end
%%%%%%%%


target=seq2;

%scale the starting point to zero
target_start=target(1,1);
%%%SA%%%
% target(:,1)=target(:,1)-target_start;
target_start_index = find(target_start == total_target(:,1));
%%%%%%%%

% scale the input such that it has the same scale as target
input_scaled = input;
input_scaled(:,1)=target(1,1)+(input(:,1)-input(1,1))*(target(end,1)-target(1,1))/...
    (input(end,1)-input(1,1));



%SA% Instead of choosing discrete ratios 1:4 - 4:1, choose ratio
%"ratio_r" from continuous data
tt= load('sedrate_dist_evenbins.txt');
marginal_dis = load('prob_original');
group_count = zeros(3,3);
group_count(1,1)=114;
group_count(1,2)=21;
group_count(1,3)=27;
group_count(2,1)=14;
group_count(2,2)=23;
group_count(2,3)=25;
group_count(3,1)=29;
group_count(3,2)=18;
group_count(3,3)=184;

%record the ratio values for each ratio
ratio_r = tt(:,2);
min_r = 0.25;
max_r = 4;

%ratio=[1,4;1,3;2,5;1,2;3,5;2,3;3,4;4,5;1,1;5,4;4,3;3,2;5,3;2,1;5,2;3,1;4,1];
rmax = length(ratio_r);

%SA% Compute probabilities of expansion(left_marginal), 1:1(mid_marginal), contraction(right_marginal)
mix_std1 = sqrt(0.0216);
mix_std2 = sqrt(0.0929);
mix_p1 = 0.642432;
mix_p2 = 1-mix_p1;
mix_mu1 = 0.0198;
mix_mu2 = -0.0297;

%0.642432*normpdf(q,0.0198, sqrt(0.0216))+0.357568*normpdf(q,-0.0297,sqrt(0.0929))
%  density_mixture_gaussian(x,mix_mu1,mix_std1,mix_mu2,mix_std2,mix_p1,mix_p2)

tt = linspace(log10(0.25),log10(0.9220),100);
s=0;
hh=tt(2)-tt(1);
for j=1:length(tt)-1
    q=(tt(j)+tt(j+1))/2;
    %kk1=0.642432*normpdf(q,0.0198, sqrt(0.0216))+0.357568*normpdf(q,-0.0297,sqrt(0.0929));
    kk1=density_mixture_gaussian(q,mix_mu1,mix_std1,mix_mu2,mix_std2,mix_p1,mix_p2);
    s=s+abs(kk1);
end
left_marginal=s*hh;

tt = linspace(log10(1.0850),log10(4),100);
s=0;
hh=tt(2)-tt(1);
for j=1:length(tt)-1
    q=(tt(j)+tt(j+1))/2;
    %kk1=0.642432*normpdf(q,0.0198, sqrt(0.0216))+0.357568*normpdf(q,-0.0297,sqrt(0.0929));
    kk1=density_mixture_gaussian(q,mix_mu1,mix_std1,mix_mu2,mix_std2,mix_p1,mix_p2);
    s=s+abs(kk1);
end
right_marginal=s*hh;

tt = linspace(log10(0.9220),log10(1.0850),100);
s=0;
hh=tt(2)-tt(1);
for j=1:length(tt)-1
    q=(tt(j)+tt(j+1))/2;
    %kk1=0.642432*normpdf(q,0.0198, sqrt(0.0216))+0.357568*normpdf(q,-0.0297,sqrt(0.0929));
    kk1=density_mixture_gaussian(q,mix_mu1,mix_std1,mix_mu2,mix_std2,mix_p1,mix_p2);
    s=s+abs(kk1);
end
mid_marginal=s*hh;

log_left_marginal = log(left_marginal);
log_right_marginal = log(right_marginal);
log_mid_marginal = log(mid_marginal);


%%%SA%%% Initial mean and standard deviation of emission
% sigma and mu are defined when we define bin_mean.
% sigma = 0.25;
% mu = 0;
%%%%%%%%


%SA% pi_nomatch remains 0
pi_nomatch = 0;


%SA% The following four variables are not used any longer.
n1=size(input_scaled,1);
n2=size(target,1);
beginpoint=target(1);
endpoint=target(n2);


%SA% Construct a discretized target such that the length of the
%discretized target's interval is 1/4 of 10th percentile of the lenth
%of input's interval.
%hh is the minimum distance between input points
hh=inf;
dist = zeros(length(input_scaled)-1,1);
for i=2:length(input_scaled)
    dist(i-1)=input_scaled(i,1)-input_scaled(i-1,1);
    if input_scaled(i,1)-input_scaled(i-1,1)<hh
        hh=input_scaled(i,1)-input_scaled(i-1,1);
    end
end

[kk1, kk2] = sort(dist);

%select bin width sunch that each bin contains at least one data point
h3 = target(2,1)-target(1,1);
for i=3:length(target)
    if target(i,1)-target(i-1,1)>h3
        h3 = target(i,1)-target(i-1,1);
    end
end

h = kk1(max(floor(length(kk1)*0.1),1))/4;
%%%SA%%%
% Set num_tar to be smaller than (target(end,1)-target(1,1))*10  
num_tar = min(target(end,1)/h, (target(end,1)-target(1,1))*10);
% num_tar = target(end,1)/h;
%%%%%%%%

%h = target(end,1)/ceil(num_tar);

% discrete_target = time axis of a discretized target
discrete_target = linspace(target(1,1),target(end,1),(ceil(num_tar)+1)/qqqq)';
h = discrete_target(2)-discrete_target(1);

%SA% Choose the largest permitted gap for an input (gap_con) and a target
% (gap_con_tar)
gap_con = ceil(0.25*length(input));
%%%SA%%%
gap_con_tar = ceil(0.25*length(discrete_target));
% gap_con_tar = ceil(10*h3/h);
%%%%%%%%

% discrete_target=[];
% discrete_target(1,1)=target(1,1);
% ss=1;
% while discrete_target(ss,1)<target(end,1)
%    ss=ss+1;
%    discrete_target(ss,1)=discrete_target(ss-1,1)+h;
% end


%%%SA%%%
% bin_mean = delO18 of a discretized target, computed by linear
% interpolation

% Estimated by step funcionts
% bin_mean_index remembers which target is used to estimate bin_mean
% For the profile HMM, we define mean and std for each target point
% (total_mu, total_sigma). For each input mu and sigma are defined in
% the same way with bin_mean
bin_mean = zeros(length(discrete_target),1);
bin_mean_index = zeros(length(discrete_target), 1);
mu = zeros(length(discrete_target),1);
sigma = zeros(length(discrete_target),1);

target_mid = target(:,1) + [diff(target(:,1))/2; 0];
tmp = find(discrete_target <= target_mid(1));
k = target_start_index;
% bin_mean(tmp) = target(1,2);
bin_mean_index(tmp) = k;
% mu(tmp) = total_mu(k);
% sigma(tmp) = total_sigma(k);
for i = 2 : length(target)
    tmp = find(discrete_target > target_mid(i-1) & discrete_target <= target_mid(i));
    
    k = i + target_start_index - 1;
    %     bin_mean(tmp) = target(i,2);
    bin_mean_index(tmp) = k;
    %     mu(tmp) = total_mu(k);
    %     sigma(tmp) = total_sigma(k);
end
% mean_mu = mean(mu);
% mean_sigma = mean(sigma);

if strcmp(target_est, 'step')
    target_mid = target(:,1) + [diff(target(:,1))/2; 0];
    tmp = find(discrete_target <= target_mid(1));
    k = target_start_index;
    bin_mean(tmp) = target(1,2);
    bin_mean_index(tmp) = k;
    mu(tmp) = total_mu(k);
    sigma(tmp) = total_sigma(k);
    for i = 2 : length(target)
        tmp = find(discrete_target > target_mid(i-1) & discrete_target <= target_mid(i));
        
        k = i + target_start_index - 1;
        bin_mean(tmp) = target(i,2);
        bin_mean_index(tmp) = k;
        mu(tmp) = total_mu(k);
        sigma(tmp) = total_sigma(k);
    end
    mean_mu = mean(mu);
    mean_sigma = mean(sigma);
end

mean_input = mean(input_scaled(:,2));
meanShift = mean_mu - mean_input;
mu = mu - meanShift;
mean_mu = mean_mu - meanShift;

% Estimation by linear interpolations
% bin_mean_index remembers which targets are used to interpolate
% bin_mean and the ratio of interpolation
% For the profile HMM, we define mean and std for each target point
% (total_mu, total_sigma). For each input mu and sigma are defined in
% the same way with bin_mean
if strcmp(target_est, 'linear')
    bin_mean = zeros(length(discrete_target),1);
    bin_mean_index_p = zeros(length(discrete_target),2);
    bin_mean_index_n = zeros(length(discrete_target),2);
    mu = zeros(length(discrete_target),1);
    sigma = zeros(length(discrete_target),1);
    
    bin_mean(1)=target(1,2);
    bin_mean_index_p(1,:) = [target_start_index 1];
    bin_mean_index_n(1,:) = [target_start_index+1 0];
    mu(1) = total_mu(target_start_index);
    sigma(1) = total_sigma(target_start_index);
    
    for i=2:length(target)
        tmp = find(discrete_target>target(i-1,1) & discrete_target<=target(i,1));
        
        k = i + target_start_index - 1;
        bin_mean_index_p(tmp,1) = k-1;
        bin_mean_index_n(tmp,1) = k;
        
        lambda = (target(i,1)-discrete_target(tmp))/(target(i,1)-target(i-1,1));
        bin_mean_index_p(tmp,2) = lambda;
        bin_mean_index_n(tmp,2) = 1 - lambda;
        
        bin_mean(tmp) = lambda*target(i-1,2) + (1-lambda)*target(i,2);
        mu(tmp) = lambda*total_mu(k-1) + (1-lambda)*total_mu(k);
        sigma(tmp) = lambda.*lambda*total_sigma(k-1) + (1-lambda).*(1-lambda)*total_sigma(k);
    end
    mean_mu = mean(mu);
    mean_sigma = mean(sigma);
end

% HMM-matching
% bin_mean = zeros(length(discrete_target),1);
% bin_mean(1)=target(1,2);
% for i=2:length(target)
%     tmp = find(discrete_target>target(i-1,1) & discrete_target<=target(i,1));
%     for j=1:length(tmp)
%         bin_mean(tmp(j))=target(i-1,2)+(target(i,2)-target(i-1,2))*...
%             (discrete_target(tmp(j))-target(i-1,1))/(target(i,1)-target(i-1,1));
%     end
% end
%%%%%%%%


% %discretize the target into corresponding bins
% discrete_target = (target(:,1):h:target(end,1))';
% interval_target_size=zeros(length(discrete_target),1);  %number of points in each bin
% interval_target = cell(length(discrete_target),1);  %record the index for points in each bin
% a=0;
% for i=1:length(interval_target_size)
%     temp=target(a+1:size(target,1),1);
%     ttt=find(temp<beginpoint+i*h)+a;
%     if length(ttt)~=0
%         a=max(ttt);
%     end;
%     interval_target{i}=ttt;
%     interval_target_size(i)=length(ttt);
%     clear ttt
% end;


%assign the index of bin for input when putting input in discreted target bins
input_bin_index = zeros(length(input),1);
ss =1;
for i = 1: length(discrete_target)
    while input_scaled(ss,1)<=discrete_target(i)
        input_bin_index(ss)=i;
        ss = ss+1;
        if ss>length(input)
            break
        end
    end
end
input_bin_index(ss:end)=length(discrete_target);


%record the number of bin distance betwen input point with previous point
%set the distance for the first element as 0
dis_bin_input = zeros(length(input),1);
for i=2:length(dis_bin_input)
    dis_bin_input(i) = input_bin_index(i) - input_bin_index(i-1);
end


input_scaled_new = input_scaled;
input_bin_index_new = input_bin_index;
dis_bin_input_new = dis_bin_input;
tt = find(dis_bin_input==0);
tt = setdiff(tt,[1]);
input_track = input(:,1);
input_track_2 = input(:,2);

for i=1:length(tt)
    input_scaled_new(tt(i)-1,1:2) = (input_scaled(tt(i),1:2)+input_scaled(tt(i)-1,1:2))/2;
    input_track(tt(i)-1)  =(input(tt(i),1)+input(tt(i)-1,1))/2;
    input_track_2(tt(i)-1)=(input(tt(i),2)+input(tt(i)-1,2))/2;
end

tmp1 = input_track(:,1);
tmp1(tt)=[];
input_track = tmp1;

tmp1 = input_scaled_new(:,1);
% tmp2=input_scaled_new(:,2);
tmp1(tt)=[];
% tmp2(tt)=[];

tmp2=input_track_2(:,1);
tmp2(tt)=[];
input_track_2=tmp2;

input_bin_index_new(tt)=[];
dis_bin_input_new(tt)=[];
% for i=1:length(tt)
%    tmp1(tt(i))=[];
%    tmp2(tt(i))=[];
%    input_bin_index_new(tt(i))=[];
%    dis_bin_input_new(tt(i))=[];
% end

dis_bin_input_new = dis_bin_input_new+1;
%%%SA%%%
% Set the size of searching window of sedimentation rate
max_dis_bin_input_new = max(dis_bin_input_new);
sorted_dis_bin_input_new = sort(dis_bin_input_new);
prac_max_dis_bin_input_new = max(4*sorted_dis_bin_input_new(ceil(length(dis_bin_input_new)*0.9)), max_dis_bin_input_new);

% memory_use_log_f = length(input_scaled_new)*length(discrete_target)*4*max_dis_bin_input_new*8/1024/1024/1024;
% memory_limit = 128;
% if memory_use_log_f < memory_limit
prac_max_dis_bin_input_new = 4*max_dis_bin_input_new;
% end
% log_f = zeros(length(input_scaled_new),length(discrete_target), prac_max_dis_bin_input_new);
%%%%%%%%

input_scaled_new = zeros(length(tmp1),2);
input_scaled_new(:,1)=tmp1;
input_scaled_new(:,2)=tmp2;


% %calculate the mean value for each bin
% bin_mean = zeros(length(discrete_target),1);
% for i=1:length(discrete_target)
%     ll1 = interval_target{i}(1);
%     ll2 = interval_target{i}(end);
%     bin_mean(i) = sum(target(ll1:ll2,2))/(ll2-ll1+1);
% end


%look up table -- record the log(pdf) and residue if data point i in input
%is aligned to j^th bin in target
% log_grid = zeros(length(input_scaled_new),length(bin_mean));
delta_grid = zeros(length(input_scaled_new),length(bin_mean));
for i=1:length(input_scaled_new)
    for j=1:length(bin_mean)
        %log_grid(i,j) = log(normpdf(input_scaled_new(i,2)-bin_mean(j),mu,sigma));
        delta_grid(i,j) = input_scaled_new(i,2)-bin_mean(j);
    end
end

%%%SA%%%
% log_grid = log(normpdf(delta_grid, mu,sigma));
log_grid = log(normpdf(delta_grid, ones(length(input_scaled_new),1)*mu', ones(length(input_scaled_new),1)*sigma'));
%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%setting the starting parameters
%setting uniform distribution as starting for EM algorithm
pi_0 = ones(rmax,1)/rmax;

group_dis = zeros(3,3);
for i=1:3
    group_dis(i,:) = group_count(i,:)/sum(group_count(i,:));
end

A_prior = zeros(rmax,rmax);
g11 = group_dis(1,1)*marginal_dis(1:8)/sum(marginal_dis(1:8));
g13 = group_dis(1,3)*marginal_dis(10:17)/sum(marginal_dis(10:17));
for i=1:8
    A_prior(i,1:8) = g11;
    A_prior(i,9) = group_dis(1,2);
    A_prior(i,10:17) = g13;
end
A_prior(9,1:8) = group_dis(2,1)*marginal_dis(1:8)/sum(marginal_dis(1:8));
A_prior(9,9) = group_dis(2,2);
A_prior(9,10:17) = group_dis(2,3)*marginal_dis(10:17)/sum(marginal_dis(10:17));

g31 = group_dis(3,1)*marginal_dis(1:8)/sum(marginal_dis(1:8));
g33 = group_dis(3,3)*marginal_dis(10:17)/sum(marginal_dis(10:17));
for i=10:17
    A_prior(i,1:8) = g31;
    A_prior(i,9) = group_dis(3,2);
    A_prior(i,10:17) = g33;
end

A = A_prior;
log_A = log(A);

%transitions from begin node to x, rate, y
pi_begin = [0.3 0.4 0.3];

%transitions from x^begin , y^begin, y^end, x^end
tao_x_begin = [0.5 0.5];
tao_y_begin = [0.5 0.5];
tao_x_end = [0.5 0.5];
tao_y_end = [0.5 0.5];

log_pi_begin = log(pi_begin);       % [I M D]
log_tao_x_begin = log(tao_x_begin); % [I M]
log_tao_y_begin = log(tao_y_begin); % [D M]
log_tao_x_end = log(tao_x_end);     % [I E]
log_tao_y_end = log(tao_y_end);     % [D E]
% log_delta_r = log(delta_r);


%%%SA%%%
delta_match = [0.9 0.1];      % M -> M or ME
delta_end = [0.4 0.3 0.3];    % ME -> E, I, or D
log_delta_match = log(delta_match);
log_delta_end = log(delta_end);

% %transitions from the rate state
% delta_r = [0.9 0.05 0.05];
% log_delta_r = log(delta_r);
%%%%%%%%


%%%SA%%%
% mean_tar should be defined outside spmd for profile HMM
log_input_gap_raw = zeros(length(input_scaled_new),1);
for i=1:length(input_scaled_new)
    log_input_gap_raw(i) = log(normpdf(input_scaled_new(i,2)-mean_tar, mean_mu , mean_sigma));
end
% mean_tar=mean(target(:,2));
% log_input_gap_raw = zeros(length(input_scaled_new),1);
% for i=1:length(input_scaled_new)
%     log_input_gap_raw(i) = log(normpdf(input_scaled_new(i,2)-mean_tar,mu,sigma));
% end
%%%%%%%%

%record starting parameters as old
mu_old = mu;
sigma_old = sigma;
pi_begin_old = pi_begin;
tao_x_begin_old = tao_x_begin;
tao_x_end_old = tao_x_end;
tao_y_begin_old = tao_y_begin;
tao_y_end_old = tao_y_end;
%%%SA%%%
delta_end_old = delta_end;
delta_match_old = delta_match;
% delta_r_old = delta_r;
log_fE_old = 0;
%%%%%%%%

%%%SA%%%
log_transition_prob_track = zeros(length(input_scaled_new),prac_max_dis_bin_input_new);
transition_prob_track = zeros(length(input_scaled_new),prac_max_dis_bin_input_new);
% log_transition_prob_track = zeros(length(input_scaled_new),max_dis_bin_input_new*4);
% transition_prob_track = zeros(length(input_scaled_new),max_dis_bin_input_new*4);
%%%%%%%%
for i=2:size(transition_prob_track,1)
    for j=1:size(transition_prob_track,2)
        
        %%%SA%%%
        % rr = dis_bin_input_new(i)/(j);
        rr = (dis_bin_input_new(i)-1)/max(j-1,0.5);
        if discreteRatio == true
            if logComp == true
                [~, rTemp] = min(abs(log10(ratio_r) - log10(rr)));
            else
                [~, rTemp] = min(abs(ratio_r - rr));
            end
            rr = ratio_r(rTemp);
        end
        %%%%%%%%
        log_rr = log10(rr);
        transition_prob_track(i,j)=density_mixture_gaussian(log_rr,mix_mu1,mix_std1,mix_mu2,mix_std2,mix_p1,mix_p2);
        log_transition_prob_track(i,j)=log(transition_prob_track(i,j));
        
    end
end

log_group_transition = log(group_dis);
log_group_transition(:,1)=log_group_transition(:,1)-log_left_marginal;
log_group_transition(:,2)=log_group_transition(:,2)-log_mid_marginal;
log_group_transition(:,3)=log_group_transition(:,3)-log_right_marginal;

group_transition = exp(log_group_transition);

%%%SA%%%
ratio_track = zeros(max_dis_bin_input_new,prac_max_dis_bin_input_new);
% ratio_track = zeros(max_dis_bin_input_new,4*max_dis_bin_input_new);
for i=1:max_dis_bin_input_new
    %     for j=1:min(4*i, prac_max_dis_bin_input_new)
    for j=1:4*i
        if logComp == true
            [~, ratio_track(i,j)] = min(abs( log10(ratio_r) - log10((i-1)/max(j-1,0.5)) ));
        else
            [~, ratio_track(i,j)] = min(abs(ratio_r-(i-1)/max(j-1,0.5)));
        end
    end
end
%%%%%%%%

% ratio_track_1=zeros(max_dis_bin_input_new,1);
% for i=1:max_dis_bin_input_new
%     %if i/0.5>=min(ratio_r) && i/0.5<=max(ratio_r)
%     [~, ratio_track_1(i)]=min(abs(ratio_r-i/0.5));
%     %end
% end

%%%SA%%%
ratio_track_new = zeros(length(input_scaled_new),prac_max_dis_bin_input_new);
% ratio_track_new = zeros(length(input_scaled_new),4*max_dis_bin_input_new);
for i=2:length(input_scaled_new)
    for j=1:prac_max_dis_bin_input_new
        %     for j=1:4*max_dis_bin_input_new
        if j~=1
            ratio_track_new(i,j) = (dis_bin_input_new(i)-1)/(j-1);
        else
            ratio_track_new(i,j)=(dis_bin_input_new(i)-1)/0.5;
        end
        if discreteRatio == true
            if logComp == true
                [~, rTemp] = min(abs(log10(ratio_r) - log10(ratio_track_new(i,j)) ));
            else
                [~, rTemp] = min(abs(ratio_r - ratio_track_new(i,j)));
            end
            ratio_track_new(i,j) = ratio_r(rTemp);
        end
    end
end
%%%%%%%%

fprintf('Preprocessing is done.\t\t\t')
toc