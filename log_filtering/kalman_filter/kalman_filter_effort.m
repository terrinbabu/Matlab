clc;
clear all;
close all;
roundn = @(x,n) round(x*10^n)./10^n;

%% input position,velocity and acceleration

mjs_value = importdata('mjs_j1_v5_25_50_75_P10D10.txt');

start_trim = 0;     
end_trim = 0;

js_time = mjs_value.data(start_trim+1:end-end_trim,2);
tau = roundn(mjs_value.data(start_trim+1:end-end_trim,15),3);

eval_js_time = js_time - js_time(1);

%% time difference & acceleration

t = 0.008;
sd3_tau = 3*25; %sd_tau25 = 15, sd_tau75 = 27;

%% Kalman filter for effort

error_zero_tolerance = 0; % to avoid constant change in control variable at low error values
P_gain = 20;
I_gain = 0;
D_gain = 0.1;

A = 1;  % A and B transformation matrix from equation considering the control variable

B = t;

P_k_minus1 = sd3_tau;  % intial process covariance matrix

R = P_k_minus1; % measurement covariance matrix

X_k_minus1 = tau(1);  % intial state - effort


delta = 0; % control variable
EF=0; % error function
last_EF = 0; % previous error function

for i = 1:size(tau)
    
    X_kp = A*X_k_minus1 + B*delta; % calculating the predicted state
    P_kp = A*P_k_minus1*A.'; % predicted process covariance matrix
    K = P_kp/(P_kp+R); % calculating the Kalman Gain
    
    Y_k = tau(i); % measured (observation) state
   
    X_k = X_kp + K*(Y_k-X_kp); % calculate the current state using the kalman filter
    P_k = (eye(1)-K)*P_kp; % updating the process covariance matrix
    
  % control variable
  
   EF = Y_k - X_k; % error function is the variation in the measured and current tau
   error(i) = EF;
    
   if abs(EF) < error_zero_tolerance
           EF = 0;
   end
   
   diff_EF = (EF - last_EF)/t;
   inter_EF = ((EF - last_EF)/2)*t;
   last_EF = EF;
   delta = P_gain*EF + I_gain*inter_EF + D_gain*diff_EF;   
   

   % state and  update 
    
    X_k_minus1 = X_k;
    P_k_minus1 = P_k;

    % for plotting the results
    estimate_tau(i) = X_k;
    d(i) = delta;
end

%% ploting the results

figure('Name','measured tau');

hold on
plot(eval_js_time,tau);
plot(eval_js_time,estimate_tau);
ylabel('tau')
legend('tau','KF-tau')

suptitle('measured tau')

figure('Name','measured tau - Speed Override 5%');

time_start_limit = find(eval_js_time>75);
time_end_limit = find(eval_js_time>130);

hold on
plot(eval_js_time(time_start_limit(1):time_end_limit(1)),tau(time_start_limit(1):time_end_limit(1)));
plot(eval_js_time(time_start_limit(1):time_end_limit(1)),estimate_tau(time_start_limit(1):time_end_limit(1)));
xlabel('time')
ylabel('tau')
legend('tau','KF-tau')

suptitle('measured tau - Speed Override 5% ')


figure('Name','measured tau - Speed Override 25%');

time_start_limit = find(eval_js_time>190);
time_end_limit = find(eval_js_time>210);

hold on
plot(eval_js_time(time_start_limit(1):time_end_limit(1)),tau(time_start_limit(1):time_end_limit(1)));
plot(eval_js_time(time_start_limit(1):time_end_limit(1)),estimate_tau(time_start_limit(1):time_end_limit(1)));
xlabel('time')
ylabel('tau')
legend('tau','KF-tau')

suptitle('measured tau - Speed Override 25% ')

figure('Name','measured tau - Speed Override 50%');

time_start_limit = find(eval_js_time>256);
time_end_limit = find(eval_js_time>266);

hold on
plot(eval_js_time(time_start_limit(1):time_end_limit(1)),tau(time_start_limit(1):time_end_limit(1)));
plot(eval_js_time(time_start_limit(1):time_end_limit(1)),estimate_tau(time_start_limit(1):time_end_limit(1)));
xlabel('time')
ylabel('tau')
legend('tau','KF-tau')

suptitle('measured tau - Speed Override 50% ')

figure('Name','measured tau - Speed Override 75%');

time_start_limit = find(eval_js_time>318);
time_end_limit = find(eval_js_time>323);

hold on
plot(eval_js_time(time_start_limit(1):time_end_limit(1)),tau(time_start_limit(1):time_end_limit(1)));
plot(eval_js_time(time_start_limit(1):time_end_limit(1)),estimate_tau(time_start_limit(1):time_end_limit(1)));
xlabel('time')
ylabel('tau')
legend('tau','KF-tau')

suptitle('measured tau - Speed Override 75% ')
