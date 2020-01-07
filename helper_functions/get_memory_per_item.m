function memory_per_item = get_memory_per_item(cost,ineq_cons,N_ineq,eq_cons,N_eq)
% memory_per_item = get_memory_per_item(cost,ineq_cons,N_ineq,eq_cons,N_eq)
%
% Given matrices describing the cost and constraints of a polynomial
% optimization problem, return the number of bytes required to store a
% single item in the list that PCBA maintains when optimizing.
%
% If the POP is of dimension l, the matrix of a polynomial should be of
% size (l+1)x(N_terms) where N_terms is the number of terms in the
% polynomial. The last row stores the coefficient of each term. For
% example, the polynomial 2x_1x_2 + 3x_1^2 is represented as:
%
%   [1 1 ; % term 1 degrees
%    2 0 ; % term 2 degrees
%    2 3]  % coefficients
%
% Author: Shreyas Kousik and Bohao Zhang
% Created: 6 Jan 2020
% Updated: nope

    % get the dimension of the problem
    l = size(cost,1) - 1 ;

    % get the multi-degree of the cost
    P = max(cost(1:l,:),[],2) ;

    % get the multi-degree of the inequality constraints
    if nargin > 1
        G = max(ineq_cons(1:l,:),[],2) ;
    else
        G = zeros(l,1) ;
        N_ineq = 0 ;
    end

    % get the multi-degree of the equality constraints
    if nargin > 3
        H = max(eq_cons(1:l,:),[],2) ;
    else
        H = zeros(l,1) ;
        N_eq = 0 ;
    end

    % get the memory per item
    memory_per_item = 4 * (2*l + prod(P+1) + N_ineq*prod(G+1) + N_eq*prod(H+1)) ;
end