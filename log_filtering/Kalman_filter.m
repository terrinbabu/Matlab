clc;
clear all;
close all;

P_gain = 5.17; % 0 to 200
I_gain = 0; % 1e4 to 7e6
D_gain = 16.59; % 0 to 200

%% Dictionary

% time,position,velocity,torque,accelearation

seq = 1;
t = 2;

j1_q = 3; j1_Dq = 9;  j1_t = 15; j1_DDq = 21;
j2_q = 4; j2_Dq = 10; j2_t = 16; j2_DDq = 22;
j3_q = 5; j3_Dq = 11; j3_t = 17; j3_DDq = 23;
j4_q = 6; j4_Dq = 12; j4_t = 18; j4_DDq = 24;
j5_q = 7; j5_Dq = 13; j5_t = 19; j5_DDq = 25;
j6_q = 8; j6_Dq = 14; j6_t = 20; j6_DDq = 26;

%% input position,velocity and acceleration

mjs_value = importdata('ijt_0kg_all.txt');

q = mjs_value.data(:,j3_q);
Dq = mjs_value.data(:,j3_Dq);
DDq  = mjs_value.data(:,j3_DDq);


%% other inputs 

error_zero_tolerance = 0.000;
t = 0.008;
sd3_q = 3*0.00001;
sd3_Dq = 3*0.2;
sd3_DDq = 3*2;

%% Kalman filter for postion, velocity , acceleration and jerk

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

for i = 1:size(q)
    
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

MA_Dq = MA_filter(Dq,25);

size_estimate_Dq = size(estimate_Dq);

for i = 2:size_estimate_Dq(2) 
slope_estimate_Dq(i) = (estimate_Dq(i) - estimate_Dq(i-1))/t;
end

manual_slope = zeros(1,size_estimate_Dq(2));
manual_slope(1591:1720)=0.17*ones();
manual_slope(1920:2045)=-0.17*ones();

Dq_error_slope = abs(slope_estimate_Dq - manual_slope);
Dq_error = abs(estimate_Dq.' - Dq);

Dq_error_slope_normalize = normalize(Dq_error_slope,'range');
Dq_error_normalize = normalize(Dq_error,'range');

j = 1; k = 1;

for i = 1:size_estimate_Dq(2)
    
    if (Dq_error_slope_normalize(i) > 0.05)
    Dq_error_slope_normalize_non_zero(j) =  Dq_error_slope_normalize(i);
    j=j+1;
    end
    
    if (Dq_error_normalize(i) > 0.05)
    Dq_error_normalize_non_zero(k) =  Dq_error_normalize(i);
    k=k+1;
    end
end

Dq_error_slope_normalize_non_zero_mean = rms(Dq_error_slope_normalize_non_zero)
Dq_error_normalize_non_zero_mean = rms(Dq_error_normalize_non_zero)

obj_fun = 0.5*Dq_error_slope_normalize_non_zero_mean + 0.5*Dq_error_normalize_non_zero_mean

% 
% estimate_Dq = estimate_Dq(1:end-9);
% 
% Dq_error = abs(estimate_Dq - MA_Dq);
% 
% size_Dq_error = size(Dq_error);
% Dq_error_mean = mean(Dq_error);
% Dq_error_norm = norm(Dq_error);
% 
% j=1;
% 
% for i = 1:size_Dq_error(1) 
%     if (Dq_error(i) > 0.001)
%         Dq_error_positive(j) = Dq_error(i);
%         j=j+1;
%     end
% end
% 
% Dq_error_sum = sum(Dq_error_positive);
% Dq_error_positive_mean = mean(Dq_error_positive)
% size_Dq_error_positive = size(Dq_error_positive);
% 
% figure('Name','Dq');
% hold on
% plot(Dq_error_positive)
% plot((Dq_error_positive_mean)*ones(1, size_Dq_error_positive(2)));


%% ploting the results

figure('Name','filtered Velocity and Acceleration');

subplot(3,1,1)
hold on
plot(q);
plot(estimate_q);
xlabel('seq')
ylabel('q')
legend('q','KF-q')

subplot(3,1,2)
hold on
plot(Dq);
plot(estimate_Dq);
xlabel('seq')
ylabel('Dq')
legend('Dq','KF-Dq')

subplot(3,1,3)
hold on
plot(DDq);
plot(estimate_DDq);
xlabel('seq')
ylabel('DDq')
legend('DDq','KF-DDq')