clc;
clear all;
close all;

roundn = @(x,n) round(x*10^n)./10^n;

it_value = importdata('ijt_no_load_sinjall_v5_1.txt');
mjs_value = importdata('mjt_no_load_sinjall_v5_1.txt');

start_trim = 0;     
end_trim = 0;

it_time = it_value.data(start_trim+1:end-end_trim,2);
iDDq1  = roundn(it_value.data(start_trim+1:end-end_trim,21),3);
iDDq2  = roundn(it_value.data(start_trim+1:end-end_trim,22),3);
iDDq3  = roundn(it_value.data(start_trim+1:end-end_trim,23),3);
iDDq4  = roundn(it_value.data(start_trim+1:end-end_trim,24),3);
iDDq5  = roundn(it_value.data(start_trim+1:end-end_trim,25),3);
iDDq6  = roundn(it_value.data(start_trim+1:end-end_trim,26),3);

mjs_time = mjs_value.data(start_trim+1:end-end_trim,2);
mq1 = roundn(mjs_value.data(start_trim+1:end-end_trim,3),3);
mq2 = roundn(mjs_value.data(start_trim+1:end-end_trim,4),3);
mq3 = roundn(mjs_value.data(start_trim+1:end-end_trim,5),3);
mq4 = roundn(mjs_value.data(start_trim+1:end-end_trim,6),3);
mq5 = roundn(mjs_value.data(start_trim+1:end-end_trim,7),3);
mq6 = roundn(mjs_value.data(start_trim+1:end-end_trim,8),3);
mDq1 = roundn(mjs_value.data(start_trim+1:end-end_trim,9),3);
mDq2 = roundn(mjs_value.data(start_trim+1:end-end_trim,10),3);
mDq3 = roundn(mjs_value.data(start_trim+1:end-end_trim,11),3);
mDq4 = roundn(mjs_value.data(start_trim+1:end-end_trim,12),3);
mDq5 = roundn(mjs_value.data(start_trim+1:end-end_trim,13),3);
mDq6 = roundn(mjs_value.data(start_trim+1:end-end_trim,14),3);
mt1 = roundn(mjs_value.data(start_trim+1:end-end_trim,15),3);
mt2 = roundn(mjs_value.data(start_trim+1:end-end_trim,16),3);
mt3 = roundn(mjs_value.data(start_trim+1:end-end_trim,17),3);
mt4 = roundn(mjs_value.data(start_trim+1:end-end_trim,18),3);
mt5 = roundn(mjs_value.data(start_trim+1:end-end_trim,19),3);
mt6 = roundn(mjs_value.data(start_trim+1:end-end_trim,20),3);

eval_it_time = it_time - it_time(1);
eval_mjs_time = mjs_time - mjs_time(1);

%% input joint and time difference

mq = mq1;
mDq = mDq1;
a = iDDq1;

for i = 1:size(mq)-1
    time_diff(i) = eval_mjs_time(i+1)-eval_mjs_time(i);
end

t = mean(time_diff);

%% adding noise
% fNoise = 50;    % Frequency [Hz]
% aNoise = 0.01;  % Amplitude
% noise  = aNoise*sin(2*pi.*eval_mjs_time.*fNoise);
% mDq = mDq + noise;

%% finding SD of velocity and Acceleration

[M,I] = max(mDq);
detrend_raise_mDq = detrend(mDq(1:I));
detrend_fall_mDq = detrend(mDq(I+1:end));
detrend_mDq = [detrend_raise_mDq.',detrend_fall_mDq.'];

sd3_mDq = 3*std(detrend_mDq);
sd3_a = 3*std(a);

%% Kalman filter for postion, velocity , acceleration and jerk

jerk_estimation_delta_t = 0.05; % must be greater than 0.001
error_zero_tolerance = 0.002; % typical value = 2e-4
P_gain = 10; % typical value = 10
D_gain = 1;

A = [ 1     t  (t^2)/2
      0     1  t      
      0     0  1       ];

B = [ (t^3)/6
      (t^2)/2
      t       ];

P_k_minus1 = [0.000001 0       0
              0        sd3_mDq 0 
              0        0       sd3_a];

R = [0.000001 0       0
     0        sd3_mDq 0 
     0        0       sd3_a];

X_k_minus1 = [ mq(1)
               mDq(1)
               a(1)  ];

J = 0;
EF=0;
area_under_mDq=0;
area_under_estimate_Dq=0;
time_limit = jerk_estimation_delta_t;
j=1;k=1;
last_EF = 0;

for i = 1:size(mq)-1
    
    X_kp = A*X_k_minus1 + B*J;
    P_kp = A*P_k_minus1*A.';
    
    P_kp(1,2)=0;    P_kp(1,3)=0;
    P_kp(2,1)=0;    P_kp(2,3)=0;
    P_kp(3,1)=0;    P_kp(3,2)=0;

    K = P_kp/(P_kp+R);
    
    Y_k = [mq(i) 
           mDq(i)
           a(i)];
    
    X_k = X_kp + K*(Y_k-X_kp);
    P_k = (eye(3)-K)*P_kp; 
    
  % jerk estimation 
  
% %    if eval_mjs_time(i) < time_limit
% %         temp_time(j) = eval_mjs_time(i);
% %         temp_estimate_Dq(j) = X_k(2,1);
% %         temp_Dq(j) = Y_k(2,1);
% %         j = j+1;
% %    else 
% %        area_under_mDq = trapz(temp_time,temp_Dq);
% %        area_under_estimate_Dq = trapz(temp_time,temp_estimate_Dq);
% %        
% %        EF = area_under_mDq-area_under_estimate_Dq;
% %        
% %        if abs(EF) < error_zero_tolerance
% %            EF = 0;
% %        end
% %        
% %        clear temp_time;
% %        clear temp_estimate_Dq;
% %        clear temp_Dq;
% %        j=1;
% %        time_limit = time_limit + jerk_estimation_delta_t;
% %    end
   

   EF = Y_k(2,1) - X_k(2,1);
   
   if abs(EF) < error_zero_tolerance
           EF = 0;
   end

   diff_EF = (EF - last_EF)/t;
   last_EF = EF;
   J = P_gain*EF + D_gain*diff_EF;
   
   jerk(i) = J;
   error_fun(i) = EF;

    
   % state update 
    
    estimate_q(i)    = X_k(1,1);
    estimate_Dq(i)   = X_k(2,1);
    estimate_acc(i)  = X_k(3,1);
    
    X_k_minus1 = X_k;
    P_k_minus1 = P_k;

    
end
  
figure('Name','filtered Velocity and Acceleration');

subplot(3,1,1)
hold on
plot(eval_mjs_time,mq);
plot(eval_mjs_time(1:end-1),estimate_q);
xlabel('time')
ylabel('q')
legend('mq','Kalman filtered q')

subplot(3,1,2)
hold on
plot(eval_mjs_time,mDq);
plot(eval_mjs_time(1:end-1),estimate_Dq);
xlabel('time')
ylabel('mDq')
legend('mDq','Kalman filtered Dq')

subplot(3,1,3)
hold on
plot(eval_mjs_time(1:end-1),estimate_acc);
xlabel('time')
ylabel('acc')
legend('Kalman filtered acc')
% ylim([-0.02,0.02])
