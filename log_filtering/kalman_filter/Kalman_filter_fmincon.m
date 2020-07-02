clc;
clear all;
close all;

P_gain = 0; % 0 to 200
I_gain = 0; % 1e4 to 7e6
D_gain = 0; % 0 to 200

X0 = [P_gain,I_gain,D_gain];

% Dq_error_norm = Kalman_filter_fn(X0)

FUN = @Kalman_filter_fn;

A = []; % here only nonlinear constraints
B = [];
Aeq = []; 
Beq=[];
LB = [0 0 0];
UB = [500 7e6 200];

[X,FVAL] = fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB);

disp(['Final Objective: ' num2str(FVAL)])
disp(['P_gain = ' num2str(X(1))])
disp(['I_gain = ' num2str(X(2))])
disp(['D_gain = ' num2str(X(3))])

