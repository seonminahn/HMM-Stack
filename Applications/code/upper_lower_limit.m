%returns the 95% upper and lower bound

function [upper_95 lower_95] = upper_lower_limit(input_scaled_new, sample_ali)

%calculate upper and lower limit
for nn=1:1000
    aa_end = max(find(sample_ali(:,nn)~=0));
    if aa_end~=size(sample_ali,1)
        for i=aa_end+1:size(sample_ali,1)
            sample_ali(i,nn)=sample_ali(aa_end,nn);
        end
    end
end

upper_95 = zeros(length(input_scaled_new),1);
lower_95 = zeros(length(input_scaled_new),1);
for i=1:length(input_scaled_new)
   temp = sample_ali(i,:);
   [aa bb] = sort(temp);
   lower_95(i)=aa(25);
   upper_95(i)=aa(975);
end


aa_end = max(find(lower_95~=0));
if aa_end~=length(lower_95)
    for i=aa_end+1:length(lower_95)
        lower_95(i)=lower_95(aa_end);
    end
end
aa_end = max(find(upper_95~=0));
if aa_end~=length(upper_95)
    for i=aa_end+1:length(upper_95)
        upper_95(i)=upper_95(aa_end);
    end
end
