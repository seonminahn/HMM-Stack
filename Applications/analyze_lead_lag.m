function analyze_lead_lag(coreName1, target_depth1, coreName2, target_depth2)
%
% example:
%     analyze_lead_lag('GeoB1032_LR04age', 3, 'GeoB1214_LR04age', 2.5)
%

%% Find the results of age estimations 
coreName = coreName1;
resultNS = ['age_estimate_' coreName '_input1_iter']; 
resultS = ['../Results/' date '/age_estimate_' coreName '_input1_iter*'];
dList = dir(resultS);
for i = 1 : length(dList)
    count1(i) = sscanf(dList(i).name, [resultNS '%d']);
end
loadFN1 = [strrep(resultS, '*', num2str(max(count1))) '.mat'];

coreName = coreName2;
resultNS = ['age_estimate_' coreName '_input1_iter']; 
resultS = ['../Results/' date '/age_estimate_' coreName '_input1_iter*'];
dList = dir(resultS);
for i = 1 : length(dList)
    count2(i) = sscanf(dList(i).name, [resultNS '%d']);
end
loadFN2 = [strrep(resultS, '*', num2str(max(count2))) '.mat'];


%%
load(loadFN1)
age1 = sample_ali;
input1 = input_scaled_new(:,2);
depth1 = input_track;

load(loadFN2)
age2 = sample_ali;
input2 = input_scaled_new(:,2);
depth2 = input_track;

clearvars -except coreName1 age1 input1 target_depth1 depth1 coreName2 age2 input2 target_depth2 depth2 

%%
[~, ind1] = min(abs(depth1-target_depth1));
[~, ind2] = min(abs(depth2-target_depth2));

str1 = coreName1;
str2 = coreName2;

dataX1 = depth1;
dataX2 = depth2;


f1 = figure;
subplot(2,2,1)
plot(dataX1, input1, 'r*-')
hold on
set(gca,'ydir', 'r')
%
axis([dataX1(1) dataX1(100) 2.25 5.5])
a = axis;
tspace = (a(2)-a(1))/50;
%
title('ODP1143')
p1 = plot(dataX1(ind1), input1(ind1), 'ro');
p1.MarkerSize = 20;
t1 = text(dataX1(ind1)+tspace, input1(ind1), 'a_1');
xlabel('Depth (m)')
ylabel(['\delta^{18}O (' char(8240) ')'])


subplot(2,2,2)
plot(dataX2, input2, 'b*-')
hold on
set(gca,'ydir', 'r')
%
axis([dataX2(1) dataX2(100) 2.25 5.5])
a = axis;
tspace = (a(2)-a(1))/50;
%
title('ODP1123')
p2 = plot(dataX2(ind2), input2(ind2), 'bo');
t2 = text(dataX2(ind2)+tspace, input2(ind2), 'a_2');
p2.MarkerSize = 20;
xlabel('Depth (m)')
ylabel(['\delta^{18}O (' char(8240) ')'])

% P(X_ind1,1 > X_ind2,2 | D, S) for given ind1 and ind2 
age1_target = age1(ind1,:);
age2_target = age2(ind2,:);

count = 0;
for i = 1 : length(age2_target)
    count = count + sum(age1_target > age2_target(i));
end
oneLeadTwo = count/(length(age1_target)*length(age2_target)); 
% str_P = ['P(' str1 '_{' num2str(ind1) '} - ' str2 '_{' num2str(ind2) '}' ' > 0) = ' num2str(oneLeadTwo)];
str_P = ['P(a_1>a_2) = ' num2str(oneLeadTwo)];


% Emprical cumulative distribution of X_ind1,1 - X_ind2,2
age1_target = age1(ind1,:);
age2_target = age2(ind2,:);

age1_target_repeat = age1_target' * ones(1, length(age2_target));
age2_target_repeat = ones(length(age1_target),1) * age2_target;
age_diff = age1_target_repeat - age2_target_repeat;

sorted_age_diff = sort(age_diff(:));
age_mean = mean(sorted_age_diff);
age_median = sorted_age_diff(round(length(sorted_age_diff)*0.5));
age_975 = sorted_age_diff(round(length(sorted_age_diff)*0.975));
age_950 = sorted_age_diff(round(length(sorted_age_diff)*0.950));
age_025 = sorted_age_diff(round(length(sorted_age_diff)*0.025));
age_050 = sorted_age_diff(round(length(sorted_age_diff)*0.050));
str = {['mean: ' num2str(age_mean)] ['median: ' num2str(age_median)] ['95% interval: [' num2str(age_025) ', ' num2str(age_975) ']' ]};

% [f, x, flo, fup] = ecdf(age_diff(:));
% f2 = figure;
% subplot(1,3,1); cdfplot(age1_target); xlabel([str1 '_{' num2str(ind1) '}']); ylabel(['F(' str1 '{' num2str(ind1) '})'])
% subplot(1,3,2); cdfplot(age2_target); xlabel([str2 '_{' num2str(ind2) '}']); ylabel(['F(' str2 '{' num2str(ind2) '})'])
% 
% subplot(1,3,3); 
% [~, stats] = cdfplot(age_diff(:))
% xlabel([str1 '_{' num2str(ind1) '} - ' str2 '_{' num2str(ind2) '}'])
% ylabel(['F(' str1 '_{' num2str(ind1) '} - ' str2 '_{' num2str(ind2) '})'])
% annotation('textbox', [0.6,0.2,0.1,0.1], 'String', [str str_P], 'linestyle', 'none')
% 
% hold on
% a = axis;
% statLine = a(1):a(2);
% plot(statLine, statLine*0 + 0.5, 'k-')
% plot(statLine, statLine*0 + 0.025, 'k-')
% plot(statLine, statLine*0 + 0.975, 'k-')



%
figure(f1)
subplot(2,2,3)
h1 = histogram(age1_target);
hold on
h2 = histogram(age2_target);
xlabel('Age estimate (ka)')
ylabel('Probability')
l1 = legend(['a_1(' strrep(str1, '_', '\_') ': Depth ' num2str(dataX1(ind1)) 'm)'], ['a_2(' strrep(str2, '_', '\_') ': Depth ' num2str(dataX2(ind2)) 'm)']);
% l1 = legend(['67^{th} point of ' str1], ['54^{th} point of ' str2]);
l1.Location = 'northwest';

subplot(2,2,4)
h3 = histogram(age_diff(:));


%
h1.BinWidth = 1;
h2.BinWidth = 1;
h3.BinWidth = 1;
h1.FaceColor = 'r';
h2.FaceColor = 'b';
h1.Normalization = 'probability';
h2.Normalization = 'probability';
h3.Normalization = 'probability';

%
subplot(2,2,4)
hold on
a = axis;
hold off
statLine2 = a(3):(a(4)-a(3))/10:a(4);
p50 = plot(statLine2*0+age_median, statLine2, 'b-');
p50.LineWidth = 2;
hold on
p025 = plot(statLine2*0+age_975, statLine2, 'r-');
p975 = plot(statLine2*0+age_025, statLine2, 'r-');
p025.LineWidth = 2;
p975.LineWidth = 2;

% p95 = plot(statLine2*0+age_950, statLine2, 'r-');
% p95.LineWidth = 2;
% p05 = plot(statLine2*0+age_050, statLine2, 'r-');
% p05.LineWidth = 2;
l2 = legend('Median', '95% interval');
% l2 = legend('Median', 'One sided 95% interval');
l2.Location = 'northwest';
h3 = histogram(age_diff(:));
h3.BinWidth = 1;
h3.FaceColor = 'k';
h3.Normalization = 'probability';
ylim(a(3:4))
% xlabel(['Difference of age estimates (kyr) (' str1 ' - ' str2  ')'] )
xlabel(['Difference of age estimates (kyr) (a_1 - a_2)'] )
ylabel('Probability')
title(str_P)

l1.Box = 'off';
l2.Box = 'off';

% subplot(2,2,1); a = axis; t1 = text(a(1)-(a(2)-a(1))/20, a(3)-(a(4)-a(3))/15, 'a'); t1.FontWeight = 'bold'; t1.FontSize = 24;
% subplot(2,2,2); a = axis; t2 = text(a(1)-(a(2)-a(1))/20, a(3)-(a(4)-a(3))/15, 'b'); t2.FontWeight = 'bold'; t2.FontSize = 24;
% subplot(2,2,3); a = axis; t3 = text(a(1)-(a(2)-a(1))/20, a(4)-(a(3)-a(4))/15, 'c'); t3.FontWeight = 'bold'; t3.FontSize = 24;
% subplot(2,2,4); a = axis; t4 = text(a(1)-(a(2)-a(1))/20, a(4)-(a(3)-a(4))/15, 'd'); t4.FontWeight = 'bold'; t4.FontSize = 24;

a1 = annotation('textbox', [0.1 0.9 0.03 0.05], 'String', 'a', 'Fontsize', 24, 'FontWeight', 'bold', 'LineStyle', 'none');
a2 = annotation('textbox', [0.1 0.45 0.03 0.05], 'String', 'c', 'Fontsize', 24, 'FontWeight', 'bold', 'LineStyle', 'none');
a3 = annotation('textbox', [0.54 0.9 0.03 0.05], 'String', 'b', 'Fontsize', 24, 'FontWeight', 'bold', 'LineStyle', 'none');
a4 = annotation('textbox', [0.54 0.45 0.03 0.05], 'String', 'd', 'Fontsize', 24, 'FontWeight', 'bold', 'LineStyle', 'none');

