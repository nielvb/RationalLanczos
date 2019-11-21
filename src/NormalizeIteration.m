function [alpha,u] = normalize_iter(vhat,what)

dotp = dot(what,vhat);
u = sqrt(1/dotp);
alpha = sqrt(1/dotp)';

end
