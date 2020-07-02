clc;
clear all;
close all;

%% Dictionary

j1 = 1; 
j2 = 2; 
j3 = 3; 
j4 = 4; 
j5 = 5; 
j6 = 6;

%%
path = '/home/terrin/projects/virtual_force_sensor/files/log_filtering/zero_phase_filter/PI_trajectory/';

time = importdata(strcat(path,'PI_time.txt'));

file_PI_q = importdata(strcat(path,'PI_org_q.txt'));
file_PI_BPF_q = importdata(strcat(path,'PI_BPF_q.txt'));

file_PI_Dq = importdata(strcat(path,'PI_org_Dq.txt'));
file_PI_BPF_Dq = importdata(strcat(path,'PI_BPF_Dq.txt'));
file_PI_ZPF_Dq = importdata(strcat(path,'PI_ZPF_Dq.txt'));
file_PI_KF_Dq = importdata(strcat(path,'PI_KF_Dq.txt'));

file_PI_DDq = importdata(strcat(path,'PI_org_DDq.txt'));
file_PI_BPF_DDq = importdata(strcat(path,'PI_BPF_DDq.txt'));
file_PI_ZPF_DDq = importdata(strcat(path,'PI_ZPF_DDq.txt'));
file_PI_KF_DDq = importdata(strcat(path,'PI_KF_DDq.txt'));

file_PI_eff = importdata(strcat(path,'PI_org_eff.txt'));
file_PI_BPF_eff = importdata(strcat(path,'PI_BPF_eff.txt'));
file_PI_ZPF_eff = importdata(strcat(path,'PI_ZPF_eff.txt'));
file_PI_KF_eff = importdata(strcat(path,'PI_KF_eff.txt'));

j = j2;

PI_q = file_PI_q(:,j);
PI_BPF_q = file_PI_BPF_q(:,j);

PI_Dq = file_PI_Dq(:,j);
PI_BPF_Dq = file_PI_BPF_Dq(:,j);
PI_ZPF_Dq = file_PI_ZPF_Dq(:,j);
PI_KF_Dq = file_PI_KF_Dq(:,j);

PI_DDq = file_PI_DDq(:,j);
PI_BPF_DDq = file_PI_BPF_DDq(:,j);
PI_ZPF_DDq = file_PI_ZPF_DDq(:,j);
PI_KF_DDq = file_PI_KF_DDq(:,j);

PI_eff = file_PI_eff(:,j);
PI_BPF_eff = file_PI_BPF_eff(:,j);
PI_ZPF_eff = file_PI_ZPF_eff(:,j);
PI_KF_eff = file_PI_KF_eff(:,j);

%% PI figures
figure('Name','filtered PI q-3');

hold on
plot(time,[PI_q,PI_BPF_q]);
xlabel('time[s]')
ylabel('q')
legend('q','BPF q')

figure('Name','filtered PI Dq-3');

subplot(3,1,1)
hold on
plot(time,PI_Dq);
plot(time,PI_BPF_Dq);
xlabel('time[s]')
ylabel('Dq')
legend('Dq','BPF Dq')

subplot(3,1,2)
hold on
plot(time,PI_Dq);
plot(time,PI_KF_Dq);
xlabel('time[s]')
ylabel('Dq')
legend('Dq','KF Dq')

subplot(3,1,3)
hold on
plot(time,PI_Dq);
plot(time,PI_ZPF_Dq);
xlabel('time[s]')
ylabel('Dq')
legend('Dq','ZPF Dq')

%% PI DDq

figure('Name','filtered PI DDq-3');

subplot(3,1,1)
hold on
plot(time,PI_DDq);
plot(time,PI_BPF_DDq);
xlabel('time[s]')
ylabel('DDq')
legend('DDq','BPF DDq')

subplot(3,1,2)
hold on
plot(time,PI_DDq);
plot(time,PI_KF_DDq);
xlabel('time[s]')
ylabel('DDq')
legend('DDq','KF DDq')

subplot(3,1,3)
hold on
plot(time,PI_DDq);
plot(time,PI_ZPF_DDq);
xlabel('time[s]')
ylabel('DDq')
legend('DDq','ZPF DDq')

%% PI eff

figure('Name','filtered PI eff-3');

subplot(3,1,1)
hold on
plot(time,PI_eff);
plot(time,PI_BPF_eff);
xlabel('time[s]')
ylabel('eff')
legend('eff','BPF eff')

subplot(3,1,2)
hold on
plot(time,PI_eff);
plot(time,PI_KF_eff);
xlabel('time[s]')
ylabel('eff')
legend('eff','KF eff')

subplot(3,1,3)
hold on
plot(time,PI_eff);
plot(time,PI_ZPF_eff);
xlabel('time[s]')
ylabel('eff')
legend('eff','ZPF eff')

%% figure

figure('Name','Joint 2 Acceleration with Zero Phase Filter ');

n=10;
coeffMA = ones(1,n)/n;

PI_ZPF_DDq = filtfilt(coeffMA,1,PI_DDq);

subplot(2,1,1)
hold on
plot(time,[PI_DDq,PI_ZPF_DDq],'LineWidth',3);
xlim([0 inf])
ylim([-inf inf])
xlabel('time[s]')
ylabel('joint Acceleration [rad/s^2]')
legend('raw signal','Zero Phase filter')
set(gca,'FontSize',16)
suptitle('Joint-2 acceleration')

subplot(2,1,2)
hold on
plot(time,PI_ZPF_DDq,'r','LineWidth',2);
xlim([0 inf])
ylim([-inf inf])
xlabel('time[s]')
ylabel('joint Acceleration [rad/s^2]')
legend('Zero Phase filter')
set(gca,'FontSize',16)
suptitle('Joint-2 acceleration')




