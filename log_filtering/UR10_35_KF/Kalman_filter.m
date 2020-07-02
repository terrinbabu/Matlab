clc;
clear all;
close all;
roundn = @(x,n) round(x*10^n)./10^n;

%% input position,velocity and acceleration

mjs_value = importdata('mjs_j1_v5_25_50_75_P10D10.txt');

start_trim = 0; %4300    
end_trim = 0; %1000

js_time = mjs_value.data(start_trim+1:end-end_trim,2);
q = roundn(mjs_value.data(start_trim+1:end-end_trim,3),3);
Dq = roundn(mjs_value.data(start_trim+1:end-end_trim,9),3);
DDq  = roundn(mjs_value.data(start_trim+1:end-end_trim,21),3);

eval_js_time = js_time - js_time(1);

%% time difference & acceleration

t = 0.008;
sd3_q = 3*0.00001;
sd3_Dq = 3*0.2;
sd3_DDq = 3*2;

%% Kalman filter for postion, velocity , acceleration and jerk

error_zero_tolerance = 0.002; % to avoid constant change in jerk value at low error values
P_gain = 10;
I_gain = 0; % typical 1e5
D_gain = 10;

A = [ 1     t  (t^2)/2   
      0     1  t      
      0     0  1       ];  % A and B transformation matrix from equation of motion considering jerk

B = [ (t^3)/6
      (t^2)/2
      t       ];

P_k_minus1 = [sd3_q 0      0
              0     sd3_Dq 0 
              0     0      sd3_DDq];  % intial process covariance matrix

R = P_k_minus1; % measurement covariance matrix

X_k_minus1 = [ q(1)
               Dq(1)
               DDq(1)];  % intial state - position, velocity and acceleration


J = 0; % jerk value
EF=0; % error function
last_EF = 0; % previous error function

for i = 1:size(q)-1
    
    X_kp = A*X_k_minus1 + B*J; % calculating the predicted state
    P_kp = A*P_k_minus1*A.'; % predicted process covariance matrix
    
    P_kp(1,2)=0;    P_kp(1,3)=0; % no state dependnacy over another state
    P_kp(2,1)=0;    P_kp(2,3)=0;
    P_kp(3,1)=0;    P_kp(3,2)=0;

    K = P_kp/(P_kp+R); % calculating the Kalman Gain

    Y_k = [q(i)  
           Dq(i)
           DDq(i)]; % measured (observation) state
   
    X_k = X_kp + K*(Y_k-X_kp); % calculate the current state using the kalman filter
    P_k = (eye(3)-K)*P_kp; % updating the process covariance matrix
    
  % jerk estimation 
  
   EF = Y_k(2,1) - X_k(2,1); % error function is the variation in the measured and current velocity
   
   if abs(EF) < error_zero_tolerance
           EF = 0;
   end

   diff_EF = (EF - last_EF)/t;
   inter_EF = ((EF - last_EF)/2)*t;
   last_EF = EF;
   J = P_gain*EF + I_gain*inter_EF + D_gain*diff_EF;
   
   % state and  update 
    
    X_k_minus1 = X_k;
    P_k_minus1 = P_k;

    % for plotting the results
    estimate_q(i)    = X_k(1,1);
    estimate_Dq(i)   = X_k(2,1);
    estimate_DDq(i)  = X_k(3,1);
    
end

%% ploting the results

figure('Name','filtered Velocity and Acceleration');

subplot(3,1,1)
hold on
plot(eval_js_time,q);
plot(eval_js_time(1:end-1),estimate_q);
xlabel('time')
ylabel('q')
legend('q','KF-q')

subplot(3,1,2)
hold on
plot(eval_js_time,Dq);
plot(eval_js_time(1:end-1),estimate_Dq);
xlabel('time')
ylabel('Dq')
legend('Dq','KF-Dq')

subplot(3,1,3)
hold on
plot(eval_js_time,DDq);
plot(eval_js_time(1:end-1),estimate_DDq);
xlabel('time')
ylabel('DDq')
legend('DDq','KF-DDq')


%%   velocity plots

figure('Name','Speed Override 5%');

time_start_limit = find(eval_js_time>75);
time_end_limit = find(eval_js_time>130);

hold on
plot(eval_js_time(time_start_limit(1):time_end_limit(1)),Dq(time_start_limit(1):time_end_limit(1)));
plot(eval_js_time(time_start_limit(1):time_end_limit(1)),estimate_Dq(time_start_limit(1):time_end_limit(1)));
xlabel('time')
ylabel('Dq')
legend('Dq','KF-Dq')

suptitle('Speed Override 5% ')


figure('Name','Speed Override 25%');

time_start_limit = find(eval_js_time>190);
time_end_limit = find(eval_js_time>210);

hold on
plot(eval_js_time(time_start_limit(1):time_end_limit(1)),Dq(time_start_limit(1):time_end_limit(1)));
plot(eval_js_time(time_start_limit(1):time_end_limit(1)),estimate_Dq(time_start_limit(1):time_end_limit(1)));
xlabel('time')
ylabel('Dq')
legend('Dq','KF-Dq')

suptitle('Speed Override 25% ')

figure('Name','Speed Override 50%');

time_start_limit = find(eval_js_time>256);
time_end_limit = find(eval_js_time>266);

hold on
plot(eval_js_time(time_start_limit(1):time_end_limit(1)),Dq(time_start_limit(1):time_end_limit(1)));
plot(eval_js_time(time_start_limit(1):time_end_limit(1)),estimate_Dq(time_start_limit(1):time_end_limit(1)));
xlabel('time')
ylabel('Dq')
legend('Dq','KF-Dq')

suptitle('Speed Override 50% ')

figure('Name','Speed Override 75%');

time_start_limit = find(eval_js_time>318);
time_end_limit = find(eval_js_time>323);

hold on
plot(eval_js_time(time_start_limit(1):time_end_limit(1)),Dq(time_start_limit(1):time_end_limit(1)));
plot(eval_js_time(time_start_limit(1):time_end_limit(1)),estimate_Dq(time_start_limit(1):time_end_limit(1)));
xlabel('time')
ylabel('Dq')
legend('Dq','KF-Dq')

suptitle('Speed Override 75% ')

figure('Name','Speed Override 95%');

time_start_limit = find(eval_js_time>360.5);
time_end_limit = find(eval_js_time>362.5);

hold on
plot(eval_js_time(time_start_limit(1):time_end_limit(1)),Dq(time_start_limit(1):time_end_limit(1)));
plot(eval_js_time(time_start_limit(1):time_end_limit(1)),estimate_Dq(time_start_limit(1):time_end_limit(1)));
xlabel('time')
ylabel('Dq')
legend('Dq','KF-Dq')

suptitle('Speed Override 95% ')
    