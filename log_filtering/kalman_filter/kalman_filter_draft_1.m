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

for i = 1:size(mq)-1
    time_diff(i) = eval_mjs_time(i+1)-eval_mjs_time(i);
end

t = mean(time_diff);
 
% % figure('Name','Measured Position & Velocity');
% % 
% % subplot(2,1,1)
% % hold on
% % plot(eval_mjs_time,mq,'o');
% % xlabel('time')
% % ylabel('mq')
% % 
% % subplot(2,1,2)
% % hold on
% % plot(eval_mjs_time,mDq,'o');
% % xlabel('time')
% % ylabel('mDq')

%% adding noise

fNoise = 500;    % Frequency [Hz]
aNoise = 0.005;  % Amplitude
noise  = aNoise*sin(2*pi.*eval_mjs_time.*fNoise);
mDq = mDq + noise;

% plot(mDq)
% hold on
% plot(mDq1)
%% velocity estimation

% method 1

for i = 2:size(mq)
    
     if i < 7
         v1(i) = (mq(i)-mq(i-1))/t;
     elseif i > 6
       % c1=1;      c2=-1;      c3=0;      c4=0;       c5=0;      c6=0;      c7=0;
       % c1=3/2;    c2=-4/2;    c3=1/2;    c4=0;       c5=0;      c6=0;      c7=0;
       % c1=11/6;   c2=-18/6;   c3=9/6;    c4=-2/6;    c5=0;      c6=0;      c7=0;
       % c1=25/12;  c2=-48/12;  c3=36/12;  c4=-16/12;  c5=3/12;   c6=0;      c7=0; 
       % c1=137/60; c2=-300/60; c3=300/60; c4=-200/60; c5=75/60;  c6=-12/60; c7=0;
        c1=147/60; c2=-360/60; c3=450/60; c4=-400/60; c5=225/60; c6=-72/60; c7=10/60;   
         
         v1(i) = ((c7*mq(i-6))+(c6*mq(i-5))+(c5*mq(i-4))+(c4*mq(i-3))+(c3*mq(i-2))+(c2*mq(i-1))+(c1*mq(i)))/(t);
     end
end

 v1 = MA_filter(v1,20);

% method 2

n = 20;
v2_final = 0;
count = 0;

for i = 2:size(mq)
    
    if count < n
        v2(i) = v2_final;
        count = count + 1;
    else
        temp_v2 = (mq(i)-mq(i-n))/(n*t);
            if temp_v2 ~= 0
               v2_final = temp_v2;
               v2(i) = v2_final;
            end
        count = 0;
    end
end

% method 3

v3 = MA_filter(mDq,20);

% % figure('Name','Velocity numerical');
% % hold on
% % plot(mDq,'o')
% % plot(v1,'o')
% % plot(v2,'o')
% % plot(v3,'o')
% % legend('mDq','v1-finite','v2-step','v3-MA')

%% acceleration Estimation

% method 1

for i = 2:size(mq)
     
    a1(i) = (mDq(i)-mDq(i-1))/t;
    
%      if i < 7
%          a1(i) = (mDq(i)-mDq(i-1))/t;
%      elseif i > 6
%         c1=1;      c2=-1;      c3=0;      c4=0;       c5=0;      c6=0;      c7=0;
%        % c1=3/2;    c2=-4/2;    c3=1/2;    c4=0;       c5=0;      c6=0;      c7=0;
%        % c1=11/6;   c2=-18/6;   c3=9/6;    c4=-2/6;    c5=0;      c6=0;      c7=0;
%        % c1=25/12;  c2=-48/12;  c3=36/12;  c4=-16/12;  c5=3/12;   c6=0;      c7=0; 
%        % c1=137/60; c2=-300/60; c3=300/60; c4=-200/60; c5=75/60;  c6=-12/60; c7=0;
%        % c1=147/60; c2=-360/60; c3=450/60; c4=-400/60; c5=225/60; c6=-72/60; c7=10/60; 
%        
%          a1(i) = ((c7*mDq(i-6))+(c6*mDq(i-5))+(c5*mDq(i-4))+(c4*mDq(i-3))+(c3*mDq(i-2))+(c2*mDq(i-1))+(c1*mDq(i)))/(t);
%      end
end

% a1 = MA_filter(a1,20);

% method 2

n = 20;
a2_final = 0;
count = 0;

for i = 2:size(mq)
    
    if count < n
        a2(i) = a2_final;
        count = count + 1;
    else
        temp_a2 = (mDq(i)-mDq(i-n))/(n*t);
            if temp_a2 ~= 0
               a2_final = temp_a2;
               a2(i) = a2_final;
            end
        count = 0;
    end
end

% % % method 3 & 4
% % 
% % a3_last = 0;
% % a4_last = 0;
% % 
% % for i = 2:size(mq)
% %     
% %      u = mDq(i-1);
% %      v = mDq(i);
% %      s = mq(i)-mq(i-1);
% %      
% %      if s == 0
% %          a3(i-1) = a3_last;
% %          a4(i-1) = a4_last;
% %      else
% %         a3(i-1) = 2*(s-(u*t))/(t^2);
% %         a4(i-1) = (v^2-u^2)/(2*s);
% %      end
% %      
% %      a3_last = a3(i-1);
% %      a4_last = a4(i-1);
% % 
% % end

% % figure('Name','Acceleration numerical');
% % hold on
% % plot(a1,'o')
% % plot(a2,'o')
% % legend('a1-MA','a2-step')


%% input velocity and acceleration 

mq = mq3;
mDq = mDq3;
a = iDDq3;

%% finding the trend and de trend of velocity

[M,I] = max(mDq);
mDq_raise = mDq(1:I);
detrend_raise_mDq = detrend(mDq_raise);
trend_raise_mDq = mDq_raise - detrend_raise_mDq;

line_eq_trend_raise_mDq = polyfit(eval_mjs_time(1:I),trend_raise_mDq,1);

mDq_fall = mDq(I+1:end);
detrend_fall_mDq = detrend(mDq_fall);
trend_fall_mDq = mDq_fall - detrend_fall_mDq;

line_eq_trend_fall_mDq = polyfit(eval_mjs_time(I+1:end),trend_fall_mDq,1);

trend_mDq = [trend_raise_mDq.',trend_fall_mDq.'];
trend_mDq = trend_mDq.';
detrend_mDq = [detrend_raise_mDq.',detrend_fall_mDq.'];
detrend_mDq = detrend_mDq.';

size_mDq = size(mDq);
sd3_mDq = 3*std(detrend_mDq);
mean_detrend_mDq = mean(detrend_mDq)*ones(1, size_mDq(1));
up_3sd_detrend_mDq = (mean(detrend_mDq)+sd3_mDq)*ones(1, size_mDq(1));
down_3sd_detrend_mDq = (mean(detrend_mDq)-sd3_mDq)*ones(1, size_mDq(1));

% % % variance calculation for decimals
% % 
% % detrend_mDq1 = detrend_mDq1*1000; // make it above one
% % mean_detrend_mDq1 = mean(detrend_mDq1)*ones(1, size_mDq1(1));
% % up_var_detrend_mDq1 = [mean(detrend_mDq1)+var(detrend_mDq1)]*ones(1, size_mDq1(1));
% % down_var_detrend_mDq1 = [mean(detrend_mDq1)-var(detrend_mDq1)]*ones(1, size_mDq1(1));
% % 
% % detrend_mDq1 = detrend_mDq1/1000;
% % mean_detrend_mDq1 = mean_detrend_mDq1/1000;
% % up_var_detrend_mDq1 = up_var_detrend_mDq1/1000;
% % down_var_detrend_mDq1 = down_var_detrend_mDq1/1000;


% % figure('Name','Velocity trend');
% % 
% % hold on
% % plot(mDq);
% % plot(trend_mDq);
% % plot(detrend_mDq);
% % plot(mean_detrend_mDq);
% % plot(up_3sd_detrend_mDq);
% % plot(down_3sd_detrend_mDq);
% % 
% % xlabel('time')
% % ylabel('mDq')
% % legend('mDq','trend mDq','detrend mDq','mean detrend mDq','up 3sd detrend mDq','down 3sd detrend mDq')

%% Acceleration trend

size_a = size(a);
sd3_a = 3*std(a);
mean_a = mean(a)*ones(1, size_a(2));
up_3sd_a = (mean(a)+sd3_a)*ones(1, size_a(2));
down_3sd_a = (mean(a)-sd3_a)*ones(1, size_a(2));

%% Kalman filter for postion, velocity , acceleration and jerk

% % %jerk calculation
% % 
% % for i = 1:size(eval_mjs_time)
% %     
% %     if eval_mjs_time(i) < 4
% %         jerk(i) = 0;
% %     elseif eval_mjs_time(i) < 16
% %         jerk(i) = 0;
% %     elseif eval_mjs_time(i) < 18
% %         jerk(i) = -0.0045;
% %     elseif eval_mjs_time(i) < 30
% %         jerk(i) = 0;
% %     elseif eval_mjs_time(i) < 33
% %         jerk(i) = 0;
% %     end
% % end

jerk_estimation_delta_t = 0.05; % must be greater than 0.001
error_zero_tolerance = 0; % typical value = 2e-4
P_gain = 10; % typical value = 10
D_gain = 30;

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
  
   if eval_mjs_time(i) < time_limit
        temp_time(j) = eval_mjs_time(i);
        temp_estimate_Dq(j) = X_k(2,1);
        temp_Dq(j) = Y_k(2,1);
        j = j+1;
   else 
       area_under_mDq = trapz(temp_time,temp_Dq);
       area_under_estimate_Dq = trapz(temp_time,temp_estimate_Dq);
       
       EF = area_under_mDq-area_under_estimate_Dq;
       
       if abs(EF) < error_zero_tolerance
           EF = 0;
       end
       
       clear temp_time;
       clear temp_estimate_Dq;
       clear temp_Dq;
       j=1;
       time_limit = time_limit + jerk_estimation_delta_t;
   end
   
   diff_EF = (EF - last_EF)/t;
   J = P_gain*EF + D_gain*diff_EF;
           
   % figure data
    
    integral_mDq(i)= area_under_mDq;
    integral_estimate_mDq(i)= area_under_estimate_Dq;
    jerk(i)=J;

    error(i)=EF;
    KG_1(i) = K(1,1);
    KG_2(i) = K(2,2);
    KG_3(i) = K(3,3);
    predicted_q(i)   = X_kp(1,1);
    predicted_Dq(i)  = X_kp(2,1);
    predicted_acc(i) = X_kp(3,1);
    estimate_q(i)    = X_k(1,1);
    estimate_Dq(i)   = X_k(2,1);
    estimate_acc(i)  = X_k(3,1);
    error_in_estimate_q(i) = P_k(1,1);
    error_in_estimate_Dq(i) = P_k(2,1);
    error_in_estimate_Dq(i) = P_k(3,1);
    
    % state update 
    
    X_k_minus1 = X_k;
    P_k_minus1 = P_k;
    last_EF = EF;
    
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
ylim([-0.02,0.02])

%% Kalman filter for postion and velocity

% % A = [ 1     t
% %       0     1 ];
% %   
% % B = [ (t^2)/2
% %       t       ];
% %  
% % R = [0.000001 0
% %      0        sd3_mDq];
% % 
% % X_k_minus1 = [ mq(1)
% %                mDq(1) ];
% %                
% % P_k_minus1 = [0.1 0 
% %               0   sd3_mDq];
% % 
% % u_k = 0;
% % 
% % for i = 1:size(mq)-1
% %     
% % %     if i == 1
% % %         u_k = 0;
% % %     else
% % %         u_k = a(i);
% % %     end
% %     
% % %     if i < I
% % %         u_k = line_eq_trend_raise_mDq(1,1);
% % %     else
% % %         u_k = line_eq_trend_fall_mDq(1,1);
% % %     end     
% %     
% %     X_kp = A*X_k_minus1 + B*u_k;
% %     P_kp = A*P_k_minus1*A.';
% %     P_kp(1,2)=0;
% %     P_kp(2,1)=0;
% %     
% %     K = P_kp/(P_kp+R);
% %     
% %     Y_k = [mq(i) 
% %            mDq(i)];
% %     
% %     X_k = X_kp + K*(Y_k-X_kp);
% %     P_k = (eye(2)-K)*P_kp;
% %     
% %     KG_1(i) = K(1,1);
% %     KG_2(i) = K(2,2);
% %     predicted_q(i) = X_kp(1,1);
% %     predicted_Dq(i) = X_kp(2,1);
% %     estimate_q(i) = X_k(1,1);
% %     estimate_Dq(i) = X_k(2,1);
% %     error_in_estimate_q(i) = P_k(1,1);
% %     error_in_estimate_Dq(i) = P_k(2,1);
% %     
% %     u_k = (X_k(2,1) - X_k_minus1(2,1))/t;
% %     estimate_a(i) = u_k;
% %     
% %     X_k_minus1 = X_k;
% %     P_k_minus1 = P_k;
% % 
% % end
% %   
% % figure('Name','filtered Velocity and Acceleration');
% % 
% % subplot(2,1,1)
% % hold on
% % plot(eval_mjs_time,mDq1);
% % plot(eval_mjs_time(1:end-1),estimate_Dq);
% % plot(eval_mjs_time,v3);
% % xlabel('time')
% % ylabel('mDq')
% % legend('mDq','estimate Dq','v3-MA')
% % 
% % subplot(2,1,2)
% % hold on
% % plot(eval_mjs_time,a);
% % plot(eval_mjs_time(1:end-1),estimate_a);
% % xlabel('time')
% % ylabel('acc')
% % legend('acc','estimate acc')


% % % % %% Kalman filter for postion, velocity and acceleration
% % % % 
% % % % A = [ 1     t  (t^2)/2
% % % %       0     1  t      
% % % %       0     0  1       ];
% % % % 
% % % % 
% % % % 
% % % % P_k_minus1 = [0.000001 0       0
% % % %               0        sd3_mDq 0 
% % % %               0        0       sd3_a];
% % % % 
% % % % R = [0.000001 0       0
% % % %      0        sd3_mDq 0 
% % % %      0        0       sd3_a];
% % % % 
% % % % X_k_minus1 = [ mq(1)
% % % %                mDq(1)
% % % %                a(1)  ];
% % % %               
% % % % 
% % % %           
% % % % u_k = 0;
% % % % 
% % % % for i = 1:size(mq)-1
% % % %     
% % % %     X_kp = A*X_k_minus1;
% % % %     P_kp = A*P_k_minus1*A.';
% % % %     
% % % %     P_kp(1,2)=0;    P_kp(1,3)=0;
% % % %     P_kp(2,1)=0;    P_kp(2,3)=0;
% % % %     P_kp(3,1)=0;    P_kp(3,2)=0;
% % % % 
% % % %     K = P_kp/(P_kp+R);
% % % %     
% % % %     Y_k = [mq(i) 
% % % %            mDq(i)
% % % %            a(i)];
% % % %     
% % % %     X_k = X_kp + K*(Y_k-X_kp);
% % % %     P_k = (eye(3)-K)*P_kp;
% % % %     
% % % %     KG_1(i) = K(1,1);
% % % %     KG_2(i) = K(2,2);
% % % %     KG_2(i) = K(3,3);
% % % %     predicted_q(i)   = X_kp(1,1);
% % % %     predicted_Dq(i)  = X_kp(2,1);
% % % %     predicted_acc(i) = X_kp(3,1);
% % % %     estimate_q(i)    = X_k(1,1);
% % % %     estimate_Dq(i)   = X_k(2,1);
% % % %     estimate_acc(i)  = X_k(3,1);
% % % %     error_in_estimate_q(i) = P_k(1,1);
% % % %     error_in_estimate_Dq(i) = P_k(2,1);
% % % %     error_in_estimate_Dq(i) = P_k(3,1);
% % % %     
% % % %     u_k = (X_k(2,1) - X_k_minus1(2,1))/t;
% % % %     
% % % %     measure_estimate_a(i) = u_k;
% % % %     
% % % %     X_k_minus1 = X_k;
% % % %     P_k_minus1 = P_k;
% % % % 
% % % % end
% % % %   
% % % % figure('Name','filtered Velocity and Acceleration');
% % % % 
% % % % subplot(3,1,1)
% % % % hold on
% % % % plot(eval_mjs_time,mq1);
% % % % plot(eval_mjs_time(1:end-1),estimate_q);
% % % % xlabel('time')
% % % % ylabel('q')
% % % % legend('mq','Kalman filtered q')
% % % % 
% % % % subplot(3,1,2)
% % % % hold on
% % % % plot(eval_mjs_time,mDq1);
% % % % plot(eval_mjs_time(1:end-1),estimate_Dq);
% % % % plot(eval_mjs_time,v3);
% % % % xlabel('time')
% % % % ylabel('mDq')
% % % % legend('mDq','Kalman filtered Dq','v3-MA')
% % % % 
% % % % subplot(3,1,3)
% % % % hold on
% % % % plot(eval_mjs_time(1:end-1),estimate_acc);
% % % % plot(eval_mjs_time(1:end-1),measure_estimate_a);
% % % % xlabel('time')
% % % % ylabel('acc')
% % % % legend('Kalman filtered acc','measure numerical acc')
% % % % 
% % % % 
% % % % sd3_measure_estimate_a = 3*std(measure_estimate_a);
% % % % 
% % % % first_filter_Dq = estimate_Dq;
% % % % 
% % % % %% 2nd filter
% % % % 
% % % % mq= mq1;
% % % % mDq = mDq1;
% % % % a = measure_estimate_a;
% % % % 
% % % % A = [ 1     t
% % % %       0     1 ];
% % % %   
% % % % B = [ (t^2)/2
% % % %       t       ];
% % % %  
% % % % R = [0.000001 0
% % % %      0        sd3_mDq];
% % % % 
% % % % X_k_minus1 = [ mq(1)
% % % %                mDq(1) ];
% % % %                
% % % % P_k_minus1 = [0.1 0 
% % % %               0   sd3_mDq];
% % % % 
% % % % for i = 1:size(mq)-1
% % % %     
% % % %     if i == 1
% % % %         u_k = 0;
% % % %     else
% % % %         u_k = a(i);
% % % %     end
% % % %     
% % % %     X_kp = A*X_k_minus1 + B*u_k;
% % % %     P_kp = A*P_k_minus1*A.';
% % % %     P_kp(1,2)=0;
% % % %     P_kp(2,1)=0;
% % % %     
% % % %     K = P_kp/(P_kp+R);
% % % %     
% % % %     Y_k = [mq(i) 
% % % %            mDq(i)];
% % % %     
% % % %     X_k = X_kp + K*(Y_k-X_kp);
% % % %     P_k = (eye(2)-K)*P_kp;
% % % %     
% % % %     KG_1(i) = K(1,1);
% % % %     KG_2(i) = K(2,2);
% % % %     predicted_q(i) = X_kp(1,1);
% % % %     predicted_Dq(i) = X_kp(2,1);
% % % %     estimate_q(i) = X_k(1,1);
% % % %     estimate_Dq(i) = X_k(2,1);
% % % %     error_in_estimate_q(i) = P_k(1,1);
% % % %     error_in_estimate_Dq(i) = P_k(2,1);
% % % %     
% % % %     X_k_minus1 = X_k;
% % % %     P_k_minus1 = P_k;
% % % % 
% % % % end
% % % % 
% % % % second_filter_woMA_Dq = estimate_Dq;
% % % % 
% % % % 
% % % % 
% % % % % mq= mq1;
% % % % % mDq = MA_filter(mDq1,20);
% % % % a = MA_filter(measure_estimate_a,20);
% % % % 
% % % % A = [ 1     t
% % % %       0     1 ];
% % % %   
% % % % B = [ (t^2)/2
% % % %       t       ];
% % % %  
% % % % R = [0.000001 0
% % % %      0        sd3_mDq];
% % % % 
% % % % X_k_minus1 = [ mq(1)
% % % %                mDq(1) ];
% % % %                
% % % % P_k_minus1 = [0.1 0 
% % % %               0   sd3_mDq];
% % % % 
% % % % for i = 1:size(mq)-1
% % % %     
% % % %     if i == 1
% % % %         u_k = 0;
% % % %     else
% % % %         u_k = a(i);
% % % %     end
% % % %     
% % % %     X_kp = A*X_k_minus1 + B*u_k;
% % % %     P_kp = A*P_k_minus1*A.';
% % % %     P_kp(1,2)=0;
% % % %     P_kp(2,1)=0;
% % % %     
% % % %     K = P_kp/(P_kp+R);
% % % %     
% % % %     Y_k = [mq(i) 
% % % %            mDq(i)];
% % % %     
% % % %     X_k = X_kp + K*(Y_k-X_kp);
% % % %     P_k = (eye(2)-K)*P_kp;
% % % %     
% % % %     KG_1(i) = K(1,1);
% % % %     KG_2(i) = K(2,2);
% % % %     predicted_q(i) = X_kp(1,1);
% % % %     predicted_Dq(i) = X_kp(2,1);
% % % %     estimate_q(i) = X_k(1,1);
% % % %     estimate_Dq(i) = X_k(2,1);
% % % %     error_in_estimate_q(i) = P_k(1,1);
% % % %     error_in_estimate_Dq(i) = P_k(2,1);
% % % %     
% % % %     X_k_minus1 = X_k;
% % % %     P_k_minus1 = P_k;
% % % % 
% % % % end
% % % % 
% % % % second_filter_wMA_Dq = estimate_Dq;
% % % % 
% % % % figure('Name','filtered-2nd Position and Velocity');
% % % % 
% % % % subplot(2,1,1)
% % % % hold on
% % % % plot(eval_mjs_time,mq1);
% % % % plot(eval_mjs_time(1:end-1),estimate_q);
% % % % xlabel('time')
% % % % ylabel('q')
% % % % legend('mq','Kalman filtered q')
% % % % 
% % % % subplot(2,1,2)
% % % % hold on
% % % % plot(eval_mjs_time,mDq1);
% % % % plot(eval_mjs_time(1:end-1),estimate_Dq);
% % % % plot(eval_mjs_time,v3);
% % % % xlabel('time')
% % % % ylabel('mDq')
% % % % legend('mDq','Kalman filtered Dq','v3-MA')
% % % % 
% % % % figure('Name','All Velocities');
% % % % hold on
% % % % plot(eval_mjs_time,mDq1);
% % % % plot(eval_mjs_time,v3)
% % % % plot(eval_mjs_time(1:end-1),first_filter_Dq)
% % % % plot(eval_mjs_time(1:end-1),second_filter_woMA_Dq)
% % % % plot(eval_mjs_time(1:end-1),estimate_Dq)
% % % % xlabel('time')
% % % % ylabel('q')
% % % % legend('mq','MA Dq','first filter Dq','second filter woMA Dq','second filter wMA Dq')
% % % % 
% % % % figure('Name','All accelearation');
% % % % hold on
% % % % % plot(eval_it_time,iDDq1,'o')
% % % % plot(eval_mjs_time(1:end-1),measure_estimate_a)
% % % % plot(eval_mjs_time(1:end-1),a)
% % % % xlabel('time')
% % % % ylabel('acc')
% % % % legend('first filter a','first filter MA a')
