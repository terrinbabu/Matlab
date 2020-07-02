% optimization problem ,
 
% nonlinear objective that the optimizer attempts to minimize,
% min x1*x4(x1+x2+x3)+x3
 
% The variable values at the optimal solution are subject to (s.t.) 
% both equality and inequality constraints.subject to,

% x1*x2*x3*x4 ≥ 25 
% x1^2 + x2^2 + x3^2 + x4^2 = 40

% lower and upper bound
% 1 ≤ x1,x2,x3,x4 ≤ 5

% intial guess
% x1=1, x2=5, x3=5, x4=1


%  The fmincon function is:
%   [X,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)

% FUN - objective function defined as function - min fun(x) 
% X0 - intial guess
% A*x ≤ B - linear inequality constraints
% Aeq*x = Beq - linear equality constraints 
% LB ≤ x ≤ UB - Lower and Upper Bound 
% NONLCON - nonlinear constraints defined as a function of c(x) < 0 ceq(x) = 0

% here,
close all
clc

FUN = @(x) x(1)*x(4)*(x(1)+x(2)+x(3))+x(3);


X0 = [1,5,5,1];
A = []; % here only nonlinear constraints
B = [];
Aeq = []; 
Beq=[];
LB = 1 * ones(4);
UB = 5 * ones(4);

NONLCON = @nlcon;

[X,FVAL] = fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON);

% show final objective
disp(['Final Objective: ' num2str(FVAL)])

% print solution
disp('Solution')
disp(['x1 = ' num2str(X(1))])
disp(['x2 = ' num2str(X(2))])
disp(['x3 = ' num2str(X(3))])
disp(['x4 = ' num2str(X(4))])


function [c,ceq] = nlcon(x)
  c = 25.0 - x(1)*x(2)*x(3)*x(4);
  ceq = sum(x.^2) - 40;
end








