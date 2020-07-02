clc;
clear all;
close all;

% P_gain = 50*2*2*2*2; % 0 to 200
% I_gain = 0; % 1e4 to 7e6
% D_gain = 10*3*1.5*2*1.1; % 0 to 200

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
mjt_value = importdata('mjt_0kg_all.txt');

% 900-3500 v1a1
% 3900-4600 v2a2
% 5000-5600 v3a3
% 6000-6400 v4a4
% 6900-7300 v5a5

s = 6900;
e = 7300;
q = mjs_value.data(:,j3_q);
Dq = mjs_value.data(:,j3_Dq);
DDq  = mjs_value.data(:,j3_DDq);
tau = mjt_value.data(:,j3_t);

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

v_level = [0.01,0.02,0.1,0.11,0.3,0.31,0.5,0.51,0.8,0.81];
P = [5, 50,  100, 200, 400, 600];
D = [5, 20,  30,  45,  90,  100];

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
   
   % gain Scheduling
   
   v = abs(Y_k(2,1));
   a = abs(X_k(3,1));
   
   I_gain = 0;

%    P_gain = (314.23*(v^3)) + (-48.79*(v^2)) + (329.45*v) + 5.73;
%    D_gain = (-120.34*(v^3)) + (180.61*(v^2)) + (40.64*v) + 0.2;

    if (v<v_level(1))
        P_gain=P(1);D_gain=D(1);
        
    elseif(v<v_level(2))
        
        x  = v;
        ya = P(1);
        yb = P(2);
        xa = v_level(1);
        xb = v_level(2);
        
        P_gain=ya+(((yb-ya)/(xb-xa))*(x-xa));
        
        ya = D(1);
        yb = D(2);
        D_gain=ya+(((yb-ya)/(xb-xa))*(x-xa));
        
    elseif (v<v_level(3))
        P_gain=P(2);D_gain=D(2);
        
    elseif(v<v_level(4))
        
        x  = v;
        ya = P(2);
        yb = P(3);
        xa = v_level(3);
        xb = v_level(4);
        
        P_gain=ya+(((yb-ya)/(xb-xa))*(x-xa));
        
        ya = D(2);
        yb = D(3);
        D_gain=ya+(((yb-ya)/(xb-xa))*(x-xa));
        
    elseif (v<v_level(5))
        P_gain=P(3);D_gain=D(3);
        
    elseif(v<v_level(6))
        
        x  = v;
        ya = P(3);
        yb = P(4);
        xa = v_level(5);
        xb = v_level(6);
        
        P_gain=ya+(((yb-ya)/(xb-xa))*(x-xa));
        
        ya = D(3);
        yb = D(4);
        D_gain=ya+(((yb-ya)/(xb-xa))*(x-xa));
        
    elseif (v<v_level(7))
        P_gain=P(4);D_gain=D(4);
        
    elseif(v<v_level(8))
        
        x  = v;
        ya = P(4);
        yb = P(5);
        xa = v_level(7);
        xb = v_level(8);
        
        P_gain=ya+(((yb-ya)/(xb-xa))*(x-xa));
        
        ya = D(4);
        yb = D(5);
        D_gain=ya+(((yb-ya)/(xb-xa))*(x-xa));
        
    elseif (v<v_level(9))
        P_gain=P(5);D_gain=D(5);
        
    elseif(v<v_level(10))
        
        x  = v;
        ya = P(5);
        yb = P(6);
        xa = v_level(9);
        xb = v_level(10);
        
        P_gain=ya+(((yb-ya)/(xb-xa))*(x-xa));
        
        ya = D(5);
        yb = D(6);
        D_gain=ya+(((yb-ya)/(xb-xa))*(x-xa));
        
    else
        P_gain=P(6);D_gain=D(6);
    end
    

%    if (v<0.01)
%      P_gain=5;D_gain=5;    
%    elseif(v<0.1)
%      P_gain=50;D_gain=20;
%    elseif(v<0.3)
%      P_gain=100;D_gain=30;
%    elseif(v<0.5)
%      P_gain=200;D_gain=45;
%    elseif(v<0.8)
%      P_gain=400;D_gain=90;
%    else
%      P_gain=800;D_gain=100;
%    end
   
   J = P_gain*EF + I_gain*inter_EF + D_gain*diff_EF;
   
   % state and  update 
    
    X_k_minus1 = X_k;
    P_k_minus1 = P_k;
    
    % for plotting the results
    
   if (v<0.01 )
      X_k(2,1) = Dq(i);
      X_k(3,1) = 0;
   end
   
    estimate_q(i)    = X_k(1,1);
    estimate_Dq(i)   = X_k(2,1);
    estimate_DDq(i)  = X_k(3,1);
    
end

n=15;
coeffMA = ones(1,n)/n;
filtfilt_Dq = filtfilt(coeffMA,1,Dq);

size_filtfilt_Dq = size(filtfilt_Dq);
size_estimate_Dq = size(estimate_Dq);

for i = 2:size_filtfilt_Dq(1) 
filtfilt_DDq(i) = (filtfilt_Dq(i) - filtfilt_Dq(i-1))/t;
end

Dq_error = abs(estimate_Dq - filtfilt_Dq.');
DDq_error = abs(estimate_DDq - filtfilt_DDq);

Dq_error_normalize = normalize(Dq_error,'range');
DDq_error_normalize = normalize(DDq_error,'range');

Dq_error_weight = 0.7;
DDq_error_weight = 1 - Dq_error_weight;

weighted_sum = Dq_error_weight*Dq_error_normalize + DDq_error_weight*DDq_error_normalize;

for i = 1:size_estimate_Dq(2) 
weighted_product(i) = (DDq_error_normalize(i).^0.5) * (Dq_error_normalize(i).^0.5);
end

weighted_sum_norm = norm(weighted_sum);
weighted_product_norm = norm(weighted_product);



%% Butter filter

Fs = 125;
cut_freq=30;

Nq = Fs/2;
norm_cut_freq = cut_freq/Nq;

[b,a] = butter(8,norm_cut_freq);
BF_q = filter(b, a, q);
BF_Dq = filter(b, a, Dq);
BF_DDq = filter(b, a, DDq);
%% ploting the results

figure('Name','filtered Position');
hold on
plot(q);
plot(BF_q);
xlabel('seq')
ylabel('q')
legend('q','BF-q')

figure('Name','filtered Velocity');
hold on
plot(Dq);
plot(BF_Dq);
plot(estimate_Dq);
xlabel('seq')
ylabel('Dq')
ylabel('Dq')
legend('Dq','BF Dq','KF-Dq')

figure('Name','filtered Accelearation');
hold on
plot(DDq);
plot(BF_DDq);
plot(estimate_DDq);
xlabel('seq')
ylabel('DDq')
legend('DDq','BF-DDq','KF-DDq')

% 900-3500 v1a1
% 3900-4600 v2a2
% 5000-5600 v3a3
% 6000-6400 v4a4
% 6900-7300 v5a5

% % figure('Name','filtered Velocity');
% % subplot(2,2,1)
% % hold on
% % plot(Dq(900:3500));
% % plot(estimate_Dq(900:3500));
% % xlabel('seq')
% % ylabel('Dq')
% % ylabel('Dq')
% % legend('Dq','KF-Dq')
% % subplot(2,2,2)
% % hold on
% % plot(Dq(3900:4600));
% % plot(estimate_Dq(3900:4600));
% % xlabel('seq')
% % ylabel('Dq')
% % ylabel('Dq')
% % legend('Dq','KF-Dq')
% % subplot(2,2,3)
% % hold on
% % plot(Dq(5000:5600));
% % plot(estimate_Dq(5000:5600));
% % xlabel('seq')
% % ylabel('Dq')
% % ylabel('Dq')
% % legend('Dq','KF-Dq')
% % subplot(2,2,4)
% % hold on
% % plot(Dq(6000:6400));
% % plot(estimate_Dq(6000:6400));
% % xlabel('seq')
% % ylabel('Dq')
% % ylabel('Dq')
% % legend('Dq','KF-Dq')



v = (0:0.001:1);
size_v = size(v);
v_level = [0.01,0.02,0.1,0.11,0.3,0.31,0.5,0.51,0.8,0.81];
P = [5, 50,  100, 200, 400, 600];
D = [5, 20,  30,  45,  90,  100];

for i = 1:size_v(2)
    
    if (v(i)<v_level(1))
        line_P(i)=P(1);line_D(i)=D(1);
        
    elseif(v(i)<v_level(2))
        
        x  = v(i);
        ya = P(1);
        yb = P(2);
        xa = v_level(1);
        xb = v_level(2);
        
        line_P(i)=ya+(((yb-ya)/(xb-xa))*(x-xa));
        
        ya = D(1);
        yb = D(2);
        line_D(i)=ya+(((yb-ya)/(xb-xa))*(x-xa));
        
    elseif (v(i)<v_level(3))
        line_P(i)=P(2);line_D(i)=D(2);
        
    elseif(v(i)<v_level(4))
        
        x  = v(i);
        ya = P(2);
        yb = P(3);
        xa = v_level(3);
        xb = v_level(4);
        
        line_P(i)=ya+(((yb-ya)/(xb-xa))*(x-xa));
        
        ya = D(2);
        yb = D(3);
        line_D(i)=ya+(((yb-ya)/(xb-xa))*(x-xa));
        
    elseif (v(i)<v_level(5))
        line_P(i)=P(3);line_D(i)=D(3);
        
    elseif(v(i)<v_level(6))
        
        x  = v(i);
        ya = P(3);
        yb = P(4);
        xa = v_level(5);
        xb = v_level(6);
        
        line_P(i)=ya+(((yb-ya)/(xb-xa))*(x-xa));
        
        ya = D(3);
        yb = D(4);
        line_D(i)=ya+(((yb-ya)/(xb-xa))*(x-xa));
        
    elseif (v(i)<v_level(7))
        line_P(i)=P(4);line_D(i)=D(4);
        
    elseif(v(i)<v_level(8))
        
        x  = v(i);
        ya = P(4);
        yb = P(5);
        xa = v_level(7);
        xb = v_level(8);
        
        line_P(i)=ya+(((yb-ya)/(xb-xa))*(x-xa));
        
        ya = D(4);
        yb = D(5);
        line_D(i)=ya+(((yb-ya)/(xb-xa))*(x-xa));
        
    elseif (v(i)<v_level(9))
        line_P(i)=P(5);line_D(i)=D(5);
        
    elseif(v(i)<v_level(10))
        
        x  = v(i);
        ya = P(5);
        yb = P(6);
        xa = v_level(9);
        xb = v_level(10);
        
        line_P(i)=ya+(((yb-ya)/(xb-xa))*(x-xa));
        
        ya = D(5);
        yb = D(6);
        line_D(i)=ya+(((yb-ya)/(xb-xa))*(x-xa));
        
    else
        line_P(i)=P(6);line_D(i)=D(6);
    end
    
end


% coeff_P = polyfit(v,P,3);
% coeff_D = polyfit(v,D,3);
% 
% x=(0:0.01:1);
% size_x = size(x);
% for i = 1:size_x(2) 
% line_P(i) = (coeff_P(1)*(x(i)^3)) + (coeff_P(2)*(x(i)^2)) + (coeff_P(3)*x(i)) + coeff_P(4);
% line_D(i) = (coeff_D(1)*(x(i)^3)) + (coeff_D(2)*(x(i)^2)) + (coeff_D(3)*x(i)) + coeff_D(4);
% end


figure('Name','v Vs P and D');

yyaxis left
plot(v,line_P);
xlabel('Velocity [m/s]')
ylabel('Proportional gain values')
ylim([0,650])

yyaxis right
plot(v,line_D,'--');
ylabel('Derivative gain values')
ylim([0,120])

legend('Proportional gain','Derivative gain')
set(gca,'FontSize',22)
suptitle('Gain Scheduling')
