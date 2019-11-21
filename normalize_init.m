function [sigma,gamma_1,r,d_1] = normalize_init(vhat_2,what_2)

dotp = dot(what_2,vhat_2);
d_1 = sqrt(dotp);
gamma_1 = sqrt(dotp)';
r = sqrt((1/dotp)*(1/(gamma_1'*d_1)));
sigma = sqrt((1/dotp)*(1/(gamma_1'*d_1)))';

end
