function k = catigorical_random_log(log_p)

u = rand;
s = 0;
p = zeros(length(log_p),1);
% for tt=1:length(log_p)
%     ssum=0;
%         for t=1:length(log_p)
%             ssum=ssum+exp(log_p(t)-log_p(tt));
%             
%         end;
%         p(tt)=1/ssum;
%     
% end

tt = max(log_p);
sss=0;
for i=1:length(log_p)
   
        sss=sss+exp(log_p(i)-tt);
    
end;
sss=tt+log(sss);
p = exp(log_p-sss);


for i=1:length(p)
   s = s+p(i);
   if u<=s
       k=i;
       break
   end
end
