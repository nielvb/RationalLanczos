function [alpha,u] = NormalizeIteration(vhat,what)
%NormalizeIteration computes normalization constants
%   Normalization constants alpha, u are computed such that
%   dot(what*alpha,vhat*u) = 1 
% INPUT
%       vhat = some vector
%       what = some vector of same length as vhat
% OUTPUT
%       alpha, u = constants such that dot(what/alpha,vhat/u) = 1 
dotp = dot(what,vhat);
u = sqrt(1/dotp);
alpha = sqrt(1/dotp)';

end
