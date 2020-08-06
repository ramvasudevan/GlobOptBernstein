%------------------------------------------------------%
% DIRECT Sample Call                                   %
% Purpose: Find the Global Minimum value of            %
%          the Goldstein-Price Test Function           %
%          using the DIRECT algorithm                  %
%          subject to the constraint                   %
%          x^2 + y^2 - 1 <= 0                          %
%          USES NAS STRATEGY                           %
%                                                      %
% You will need the accompanying programs              %
% Direct.m, and gp_con.m to run this code              %
% successfully.                                        %
%                                                      %
% These codes can be found at                          %
% http://www4.ncsu.edu/~definkel/research/index.html   %
%------------------------------------------------------%
clear all;
% 1. Establish bounds for variables
bounds = [-2 2;-2 2];

% 2. Send options to Direct
%    We tell DIRECT that the globalmin = 3
%    It will stop within 0.01% of solution
options.testflag  = 1;
options.globalmin = 3;
options.showits   = 1;
options.tol       = 0.01;

% NEW Update options to use implicit constraints
options.impcons = 1;

% 2a. NEW!
% Pass Function as part of a Matlab Structure
Problem.f = 'gp_con';


% 3. Call DIRECT
[fmin,x,hist] = Direct(Problem,bounds,options);

% 4. Plot iteration statistics
plot(hist(:,2),hist(:,3))
xlabel('Fcn Evals');
ylabel('f_{min}');
title('Iteration Statistics for GP test Function');
