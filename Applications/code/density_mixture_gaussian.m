function pp = density_mixture_gaussian(x,mu1,std1,mu2,std2,p1,p2)

pp = p1*normpdf(x,mu1,std1)+p2*normpdf(x,mu2,std2);

