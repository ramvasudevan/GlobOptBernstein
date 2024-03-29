function [value,fflag] = gp_con(x)
%----------------------------------------------------------
% Goldstein-Price Test Function for Nonlinear Optimization
%
% Taken from "Towards Global Optimisation 2",edited by L.C.W. Dixon and G.P.
% Szego, North-Holland Publishing Company, 1978. ISBN 0 444 85171 2
%
% -2 <= x1 <= 2
% -2 <= x2 <= 2
% fmin = 3
% xmin = 0
%       -1
%----------------------------------------------------------

%---------------------------------------------------------%
% http://www4.ncsu.edu/~definkel/research/index.html
%---------------------------------------------------------%

con_value = x(1)*x(1) + x(2)*x(2) - 1;
if con_value > 0
    fflag = 1;
    value = 0;
else
    fflag = 0;

    x1 = x(1); x2 = x(2);
    value =(1+(x1+x2+1).^2.*(19-14.*x1+3.*x1.^2-14.*x2+6.*x1.*x2+3.*x2.^2))...
     .*(30+(2.*x1-3.*x2).^2.*(18-32.*x1+12.*x1.^2+48.*x2-36.*x1.*x2+27.*x2.^2));
end
return;
