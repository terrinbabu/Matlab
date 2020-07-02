clc;
clear all;
close all;

roundn = @(x,n) round(x*10^n)./10^n;

mjs_value = importdata('mjs_j1_v5_25_50_75_KF2.txt');
fjs_value = importdata('fjs_j1_v5_25_50_75_KF2.txt');

start_trim = 0;     
end_trim = 0;

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
mDDq1  = roundn(mjs_value.data(start_trim+1:end-end_trim,21),3);
mDDq2  = roundn(mjs_value.data(start_trim+1:end-end_trim,22),3);
mDDq3  = roundn(mjs_value.data(start_trim+1:end-end_trim,23),3);
mDDq4  = roundn(mjs_value.data(start_trim+1:end-end_trim,24),3);
mDDq5  = roundn(mjs_value.data(start_trim+1:end-end_trim,25),3);
mDDq6  = roundn(mjs_value.data(start_trim+1:end-end_trim,26),3);

fjs_time = fjs_value.data(start_trim+1:end-end_trim,2);
fq1 = roundn(fjs_value.data(start_trim+1:end-end_trim,3),3);
fq2 = roundn(fjs_value.data(start_trim+1:end-end_trim,4),3);
fq3 = roundn(fjs_value.data(start_trim+1:end-end_trim,5),3);
fq4 = roundn(fjs_value.data(start_trim+1:end-end_trim,6),3);
fq5 = roundn(fjs_value.data(start_trim+1:end-end_trim,7),3);
fq6 = roundn(fjs_value.data(start_trim+1:end-end_trim,8),3);
fDq1 = roundn(fjs_value.data(start_trim+1:end-end_trim,9),3);
fDq2 = roundn(fjs_value.data(start_trim+1:end-end_trim,10),3);
fDq3 = roundn(fjs_value.data(start_trim+1:end-end_trim,11),3);
fDq4 = roundn(fjs_value.data(start_trim+1:end-end_trim,12),3);
fDq5 = roundn(fjs_value.data(start_trim+1:end-end_trim,13),3);
fDq6 = roundn(fjs_value.data(start_trim+1:end-end_trim,14),3);
ft1 = roundn(fjs_value.data(start_trim+1:end-end_trim,15),3);
ft2 = roundn(fjs_value.data(start_trim+1:end-end_trim,16),3);
ft3 = roundn(fjs_value.data(start_trim+1:end-end_trim,17),3);
ft4 = roundn(fjs_value.data(start_trim+1:end-end_trim,18),3);
ft5 = roundn(fjs_value.data(start_trim+1:end-end_trim,19),3);
ft6 = roundn(fjs_value.data(start_trim+1:end-end_trim,20),3);
fDDq1  = roundn(fjs_value.data(start_trim+1:end-end_trim,21),3);
fDDq2  = roundn(fjs_value.data(start_trim+1:end-end_trim,22),3);
fDDq3  = roundn(fjs_value.data(start_trim+1:end-end_trim,23),3);
fDDq4  = roundn(fjs_value.data(start_trim+1:end-end_trim,24),3);
fDDq5  = roundn(fjs_value.data(start_trim+1:end-end_trim,25),3);
fDDq6  = roundn(fjs_value.data(start_trim+1:end-end_trim,26),3);

eval_mjs_time = mjs_time - mjs_time(1);
eval_fjs_time = fjs_time - fjs_time(1);

%% Position 

figure ('Name','Position with Filter');

q_ylim_l = -4;
q_ylim_u = 4;

for i = 1:6
    
    subplot(3,2,i)
    plot(eval_mjs_time,roundn(mjs_value.data(start_trim+1:end-end_trim,i+2),3))
    hold on
    plot(eval_fjs_time,roundn(fjs_value.data(start_trim+1:end-end_trim,i+2),3))    
    
    xlabel('time')
    ylabel(strcat('q',int2str(i),'/fq',int2str(i)))
    legend(strcat('q',int2str(i)),strcat('fq',int2str(i)))
    ylim([q_ylim_l,q_ylim_u])
end

suptitle('Joint Postion with Filter')

%% Velocity

figure ('Name','Velocity  with Filter')

Dq_ylim_l = -0.15;
Dq_ylim_u = 0.15;

for i = 1:6
    
    subplot(3,2,i)
    plot(eval_mjs_time,roundn(mjs_value.data(start_trim+1:end-end_trim,i+8),3))
    hold on
    plot(eval_fjs_time,roundn(fjs_value.data(start_trim+1:end-end_trim,i+8),3))
    
    xlabel('time')
    ylabel(strcat('Dq',int2str(i),'/fDq',int2str(i)))
    legend(strcat('Dq',int2str(i)),strcat('fDq',int2str(i)))
%     ylim([Dq_ylim_l,Dq_ylim_u])
end

suptitle('Joint Velocity  with Filter')


figure ('Name','Joint 1 Velocity  with Filter')

plot(eval_mjs_time,roundn(mjs_value.data(start_trim+1:end-end_trim,1+8),3))
hold on
plot(eval_fjs_time,roundn(fjs_value.data(start_trim+1:end-end_trim,1+8),3))
    
xlabel('time')
ylabel('Dq1/fDq1')
legend('Dq1','fDq1')

suptitle('Joint 1 Velocity  with Filter')

%% Tau

figure ('Name','Measured tau with Filter')

for i = 1:6
    
    subplot(3,2,i)
    plot(eval_mjs_time,roundn(mjs_value.data(start_trim+1:end-end_trim,i+14),3))
    hold on
    plot(eval_fjs_time,roundn(fjs_value.data(start_trim+1:end-end_trim,i+14),3))
    
    xlabel('time')
    ylabel(strcat('tau',int2str(i),'/ftau',int2str(i)))
    legend(strcat('tau',int2str(i)),strcat('ftau',int2str(i)))
end

suptitle('Measured tau  with Filter')


figure ('Name','Joint 1 Measured tau  with Filter')

plot(eval_mjs_time,roundn(mjs_value.data(start_trim+1:end-end_trim,1+14),3))
hold on
plot(eval_fjs_time,roundn(fjs_value.data(start_trim+1:end-end_trim,1+14),3))
    
xlabel('time')
ylabel('tau1/ftau1')
legend('tau1','ftau1')

suptitle('Joint 1 Measured tau with Filter')

%% Acceleration

figure ('Name','Acceleration with Filter')

for i = 1:6
    
    
    
    subplot(3,2,i)
    plot(eval_mjs_time,roundn(mjs_value.data(start_trim+1:end-end_trim,i+20),3))
    hold on
    plot(eval_fjs_time,roundn(fjs_value.data(start_trim+1:end-end_trim,i+20),3))
    
    xlabel('time')
    ylabel(strcat('DDq',int2str(i),'/fDDq',int2str(i)))
    legend(strcat('DDq',int2str(i)),strcat('fDDq',int2str(i)))
    
end

suptitle('Joint Acceleration with Filter')

figure ('Name','Joint 1 Acceleration  with Filter')

plot(eval_mjs_time,roundn(mjs_value.data(start_trim+1:end-end_trim,1+20),3))
hold on
plot(eval_fjs_time,roundn(fjs_value.data(start_trim+1:end-end_trim,1+20),3))
    
xlabel('time')
ylabel('DDq1/fDDq1')
legend('DDq1','fDDq1')

suptitle('Joint 1 Acceleration  with Filter')
