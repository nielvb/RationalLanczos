function [sigma,gamma_1,r,d_1] = NormalizeInit(vhat,what)
%NormalizeInit  computes normalization constants
%   Normalization constants sigma, gamma_1, r, d_1 are computed such that
%   dot(what*(alpha*gamma),vhat*(r*d_1)) = 1 
% INPUT
%       vhat = some vector
%       what = some vector of same length as vhat
% OUTPUT
%       sigma, gamma_1, r, d_1 = constants such that dot(what*(alpha*gamma),vhat*(r*d_1)) = 1 
dotp = dot(what,vhat);
d_1 = sqrt(dotp);
gamma_1 = sqrt(dotp)';
r = sqrt((1/dotp)*(1/(gamma_1'*d_1)));
sigma = sqrt((1/dotp)*(1/(gamma_1'*d_1)))';

end
