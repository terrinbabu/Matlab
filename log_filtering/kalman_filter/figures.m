clc;
clear all;
close all;

%% Dictionary

% experiment

org = 1;
jskf = 2;
pikf = 3;
pikf_jskf = 4;

% ros topics

it  = 1;
ft  = 2;
et  = 3;
mjs = 4;
rt  = 5;
vt  = 6;

% time,position,velocity,torque,accelearation

seq = 1;
t = 2;

j1_q = 3; j1_Dq = 9;  j1_t = 15; j1_DDq = 21;
j2_q = 4; j2_Dq = 10; j2_t = 16; j2_DDq = 22;
j3_q = 5; j3_Dq = 11; j3_t = 17; j3_DDq = 23;
j4_q = 6; j4_Dq = 12; j4_t = 18; j4_DDq = 24;
j5_q = 7; j5_Dq = 13; j5_t = 19; j5_DDq = 25;
j6_q = 8; j6_Dq = 14; j6_t = 20; j6_DDq = 26;

%% storing all the topic values 

value = {}; % stores all the topics from all the txt files

value{org,it} = importdata('ijt_0kg_all.txt');
value{org,ft} = importdata('fjt_0kg_all.txt');
value{org,et} = importdata('ejt_0kg_all.txt');
value{org,mjs} = importdata('mjt_0kg_all.txt');
value{org,rt} = importdata('rjt_0kg_all.txt');
value{org,vt} = importdata('vf_0kg_all.txt');

value{jskf,it} = importdata('ijt_0kg_all_jskf2.txt');
value{jskf,ft} = importdata('fjt_0kg_all_jskf2.txt');
value{jskf,et} = importdata('ejt_0kg_all_jskf2.txt');
value{jskf,mjs} = importdata('mjt_0kg_all_jskf2.txt');
value{jskf,rt} = importdata('rjt_0kg_all_jskf2.txt');
value{jskf,vt} = importdata('vf_0kg_all_jskf2.txt');

value{pikf,it} = importdata('ijt_0kg_all_pikf2.txt');
value{pikf,ft} = importdata('fjt_0kg_all_pikf2.txt');
value{pikf,et} = importdata('ejt_0kg_all_pikf2.txt');
value{pikf,mjs} = importdata('mjt_0kg_all_pikf2.txt');
value{pikf,rt} = importdata('rjt_0kg_all_pikf2.txt');
value{pikf,vt} = importdata('vf_0kg_all_pikf2.txt');

value{pikf_jskf,it} = importdata('ijt_0kg_all_pikf2_jskf2.txt');
value{pikf_jskf,ft} = importdata('fjt_0kg_all_pikf2_jskf2.txt');
value{pikf_jskf,et} = importdata('ejt_0kg_all_pikf2_jskf2.txt');
value{pikf_jskf,mjs} = importdata('mjt_0kg_all_pikf2_jskf2.txt');
value{pikf_jskf,rt} = importdata('rjt_0kg_all_pikf2_jskf2.txt');
value{pikf_jskf,vt} = importdata('vf_0kg_all_pikf2_jskf2.txt');

%% sequence

% v1a1_s = 900; v1a1_e = 2300;
% v2a2_s = 2800; v2a2_e = 3100;
% v1a2_s = 3600; v1a2_e = 5000;
% v2a1_s = 5400; v2a1_e = 6300;

%% time

v1a1_s = 7.5; v1a1_e = 18.3;
v2a2_s = 22.5; v2a2_e = 25;
v1a2_s = 29.4; v1a2_e = 39.5;
v2a1_s = 43.6; v2a1_e = 50.4;

%% changing the start time stamp value to O

bags = {value{org,:},value{jskf,:},value{pikf,:},value{pikf_jskf,:}}; % same rosbag file

size_bags = size(bags);
start_time = zeros(1,size_bags(1));

for i=1:size_bags(1)
    
start_time(i) = bags{i,1}.data(1,t); % 1st time stamp 

    for j=1:size_bags(2)
        if( bags{i,j}.data(1,t)  < start_time(i) )
            start_time(i) = bags{i,j}.data(1,t); 
        end
    end
end


for i=1:size_bags(1)
    for j=1:size_bags(2)
        bags{i,j}.data(:,2) = bags{i,j}.data(:,2)-start_time(i);
    end
end

size_values = size(value);

for i=1:size_values(2)
    
value{org,i}.data(:,2) = value{org,i}.data(:,2) - start_time(1);
value{jskf,i}.data(:,2) = value{jskf,i}.data(:,2) - start_time(1);
value{pikf,i}.data(:,2) = value{pikf,i}.data(:,2) - start_time(1);
value{pikf_jskf,i}.data(:,2) = value{pikf_jskf,i}.data(:,2) - start_time(1);

end

%% Butter filter

Fs = 125;
cut_freq=30;

Nq = Fs/2;
norm_cut_freq = cut_freq/Nq;

[b,a] = butter(8,norm_cut_freq);
BF_q = filter(b, a, value{org,mjs}.data(:,j3_q));
BF_Dq = filter(b, a, value{org,mjs}.data(:,j3_Dq));
BF_DDq = filter(b, a, value{org,mjs}.data(:,j3_Dq));
BF_t = filter(b, a, value{org,mjs}.data(:,j3_t));


%% measured position,Velocity, accerleration and torque

figure ('Name','position joint 3')
hold on
plot(value{org,mjs}.data(:,t),value{org,mjs}.data(:,j3_q))
plot(value{jskf,mjs}.data(:,t),value{jskf,mjs}.data(:,j3_q))
plot(value{org,mjs}.data(:,t),BF_q)
xlabel('time [s]')
ylabel('joint position [rad]')
legend('raw signal','kalman filter','Butterworth filter')
set(gca,'FontSize',16)
suptitle('Joint-3 Postion')


figure ('Name','Velocity joint 3')
hold on
plot(value{org,mjs}.data(:,t),value{org,mjs}.data(:,j3_Dq))
plot(value{jskf,mjs}.data(:,t),value{jskf,mjs}.data(:,j3_Dq))
plot(value{org,mjs}.data(:,t),BF_Dq)
xlabel('time [s]')
ylabel('joint velocity [rad/s]')
legend('raw signal','kalman filter','Butterworth filter')
set(gca,'FontSize',16)
suptitle('Joint-3 Velocity')

figure ('Name','accerleration joint 3')
hold on
plot(value{org,it}.data(:,t),value{org,it}.data(:,j3_DDq))
plot(value{jskf,it}.data(:,t),value{jskf,it}.data(:,j3_DDq))
plot(value{org,mjs}.data(:,t),BF_DDq)
xlabel('time [s]')
ylabel('joint accerleration [rad/s^2]')
legend('raw signal','kalman filter','Butterworth filter')
set(gca,'FontSize',16)
suptitle('Joint-3 Accerleration')

figure ('Name','measured torque joint 3')
hold on
plot(value{org,mjs}.data(:,t),value{org,mjs}.data(:,j3_t))
plot(value{jskf,mjs}.data(:,t),value{jskf,mjs}.data(:,j3_t))
plot(value{org,mjs}.data(:,t),BF_t)
xlabel('time [s]')
ylabel('measured torque[Nm]')
legend('raw signal','kalman filter','Butterworth filter')
set(gca,'FontSize',16)
suptitle('Joint-3 measured torque')


%% 

size_org = size(value{org,it}.data(:,t));
count_1 = 0;
count_2 = 0;
count_3 = 0;
count_4 = 0;
count_5 = 0;
count_6 = 0;
count_7 = 0;
count_8 = 0;

for i = 1:size_org(1)
    
    if count_1 == 0
    if value{org,it}.data(i,t) > v1a1_s
        v1a1_os = i;
        count_1 = 1;
    end
    end
 

    
    if count_2 == 0
    if value{org,it}.data(i,t) > v1a1_e
        v1a1_oe = i;
        count_2 = 1;
    end
    end
    
    
    if count_3 == 0
    if value{org,it}.data(i,t) > v2a2_s
        v2a2_os = i;
        count_3 = 1;
    end
    end

    
    if count_4 == 0
    if value{org,it}.data(i,t) > v2a2_e
        v2a2_oe = i;
        count_4 = 1;
    end
    end
    
    if count_5 == 0
    if value{org,it}.data(i,t) > v1a2_s
        v1a2_os = i;
        count_5 = 1;
    end
    end

    
    if count_6 == 0
    if value{org,it}.data(i,t) > v1a2_e
        v1a2_oe = i;
        count_6 = 1;
    end
    end
    
     if count_7 == 0
    if value{org,it}.data(i,t) > v2a1_s
        v2a1_os = i;
        count_7 = 1;
    end
    end

    
    if count_8 == 0
    if value{org,it}.data(i,t) > v2a1_e
        v2a1_oe = i;
        count_8 = 1;
    end
    end
    
    
end

size_jskf = size(value{jskf,it}.data(:,t));
count_1 = 0;
count_2 = 0;
count_3 = 0;
count_4 = 0;
count_5 = 0;
count_6 = 0;
count_7 = 0;
count_8 = 0;

for i = 1:size_jskf(1)
    
    if count_1 == 0
    if value{jskf,it}.data(i,t) > v1a1_s
        v1a1_js = i;
        count_1 = 1;
    end
    end
 

    
    if count_2 == 0
    if value{jskf,it}.data(i,t) > v1a1_e
        v1a1_je = i;
        count_2 = 1;
    end
    end
    

    
    if count_3 == 0
    if value{jskf,it}.data(i,t) > v2a2_s
        v2a2_js = i;
        count_3 = 1;
    end
    end

    
    if count_4 == 0
    if value{jskf,it}.data(i,t) > v2a2_e
        v2a2_je = i;
        count_4 = 1;
    end
    end
    
    if count_5 == 0
    if value{jskf,it}.data(i,t) > v1a2_s
        v1a2_js = i;
        count_5 = 1;
    end
    end

    
    if count_6 == 0
    if value{jskf,it}.data(i,t) > v1a2_e
        v1a2_je = i;
        count_6 = 1;
    end
    end
    
     if count_7 == 0
    if value{jskf,it}.data(i,t) > v2a1_s
        v2a1_js = i;
        count_7 = 1;
    end
    end

    
    if count_8 == 0
    if value{jskf,it}.data(i,t) > v2a1_e
        v2a1_je = i;
        count_8 = 1;
    end
    end
    
end

v1a1_org_time = (value{org,it}.data(v1a1_os:v1a1_oe,t))-(value{org,it}.data(v1a1_os,t));
v1a1_org_it = (value{org,it}.data(v1a1_os:v1a1_oe,j3_t));
v1a1_org_ft = (value{org,ft}.data(v1a1_os:v1a1_oe,j3_t));
v1a1_org_et = (value{org,et}.data(v1a1_os:v1a1_oe,j3_t));
v1a1_jskf_time = (value{jskf,it}.data(v1a1_js:v1a1_je,t))-(value{jskf,it}.data(v1a1_js,t));
v1a1_jskf_it = (value{jskf,it}.data(v1a1_js:v1a1_je,j3_t));
v1a1_jskf_ft = (value{jskf,ft}.data(v1a1_js:v1a1_je,j3_t));
v1a1_jskf_et = (value{jskf,et}.data(v1a1_js:v1a1_je,j3_t));

v2a2_org_time = (value{org,it}.data(v2a2_os:v2a2_oe,t))-(value{org,it}.data(v2a2_os,t));
v2a2_org_it = (value{org,it}.data(v2a2_os:v2a2_oe,j3_t));
v2a2_org_ft = (value{org,ft}.data(v2a2_os:v2a2_oe,j3_t));
v2a2_org_et = (value{org,et}.data(v2a2_os:v2a2_oe,j3_t));
v2a2_jskf_time = (value{jskf,it}.data(v2a2_js:v2a2_je,t))-(value{jskf,it}.data(v2a2_js,t));
v2a2_jskf_it = (value{jskf,it}.data(v2a2_js:v2a2_je,j3_t));
v2a2_jskf_ft = (value{jskf,ft}.data(v2a2_js:v2a2_je,j3_t));
v2a2_jskf_et = (value{jskf,et}.data(v2a2_js:v2a2_je,j3_t));

v1a2_org_time = (value{org,it}.data(v1a2_os:v1a2_oe,t))-(value{org,it}.data(v1a2_os,t));
v1a2_org_it = (value{org,it}.data(v1a2_os:v1a2_oe,j3_t));
v1a2_org_ft = (value{org,ft}.data(v1a2_os:v1a2_oe,j3_t));
v1a2_org_et = (value{org,et}.data(v1a2_os:v1a2_oe,j3_t));
v1a2_jskf_time = (value{jskf,it}.data(v1a2_js:v1a2_je,t))-(value{jskf,it}.data(v1a2_js,t));
v1a2_jskf_it = (value{jskf,it}.data(v1a2_js:v1a2_je,j3_t));
v1a2_jskf_ft = (value{jskf,ft}.data(v1a2_js:v1a2_je,j3_t));
v1a2_jskf_et = (value{jskf,et}.data(v1a2_js:v1a2_je,j3_t));

v2a1_org_time = (value{org,it}.data(v2a1_os:v2a1_oe,t))-(value{org,it}.data(v2a1_os,t));
v2a1_org_it = (value{org,it}.data(v2a1_os:v2a1_oe,j3_t));
v2a1_org_ft = (value{org,ft}.data(v2a1_os:v2a1_oe,j3_t));
v2a1_org_et = (value{org,et}.data(v2a1_os:v2a1_oe,j3_t));
v2a1_jskf_time = (value{jskf,it}.data(v2a1_js:v2a1_je,t))-(value{jskf,it}.data(v2a1_js,t));
v2a1_jskf_it = (value{jskf,it}.data(v2a1_js:v2a1_je,j3_t));
v2a1_jskf_ft = (value{jskf,ft}.data(v2a1_js:v2a1_je,j3_t));
v2a1_jskf_et = (value{jskf,et}.data(v2a1_js:v2a1_je,j3_t));

figure ('Name','inertial torque joint 3')

subplot(2,2,1)
hold on
plot(v1a1_org_time,v1a1_org_it)
plot(v1a1_jskf_time,v1a1_jskf_it,'--','LineWidth',2)
xlabel('time [s]')
ylabel('inertial torque[Nm]')
xlim([0 inf])
ylim([-inf inf])
title('Trajectory 1')
set(gca,'FontSize',20)

subplot(2,2,2)
hold on
plot(v2a2_org_time,v2a2_org_it)
plot(v2a2_jskf_time,v2a2_jskf_it,'--','LineWidth',2)
xlabel('time [s]')
ylabel('inertial torque[Nm]')
legend('raw signal','kalman filter')
xlim([0 inf])
ylim([-inf inf])
title('Trajectory 2')
set(gca,'FontSize',20)

subplot(2,2,3)
hold on
plot(v1a2_org_time,v1a2_org_it)
plot(v1a2_jskf_time,v1a2_jskf_it,'--','LineWidth',2)
xlabel('time [s]')
ylabel('inertial torque[Nm]')
xlim([0 inf])
ylim([-inf inf])
title('Trajectory 3')
set(gca,'FontSize',20)

subplot(2,2,4)
hold on
plot(v2a1_org_time,v2a1_org_it)
plot(v2a1_jskf_time,v2a1_jskf_it,'--','LineWidth',2)
xlabel('time [s]')
ylabel('inertial torque[Nm]')
xlim([0 inf])
ylim([-inf inf])
title('Trajectory 4')
set(gca,'FontSize',20)

suptitle('Joint-3 inertial torque')

%%

figure ('Name','friction torque joint 3')

subplot(2,2,1)
hold on
plot(v1a1_org_time,v1a1_org_ft,'LineWidth',2)
plot(v1a1_jskf_time,v1a1_jskf_ft,'--','LineWidth',2)
xlabel('time [s]')
ylabel('friction torque[Nm]')
xlim([0 inf])
title('Trajectory 1')
set(gca,'FontSize',20)

subplot(2,2,2)
hold on
plot(v2a2_org_time,v2a2_org_ft,'LineWidth',2)
plot(v2a2_jskf_time,v2a2_jskf_ft,'--','LineWidth',2)
xlabel('time [s]')
ylabel('friction torque[Nm]')
legend('raw signal','kalman filter')
xlim([0 inf])
title('Trajectory 2')
set(gca,'FontSize',20)

subplot(2,2,3)
hold on
plot(v1a2_org_time,v1a2_org_ft,'LineWidth',2)
plot(v1a2_jskf_time,v1a2_jskf_ft,'--','LineWidth',2)
xlabel('time [s]')
ylabel('friction torque[Nm]')
xlim([0 inf])
title('Trajectory 3')
set(gca,'FontSize',20)

subplot(2,2,4)
hold on
plot(v2a1_org_time,v2a1_org_ft,'LineWidth',2)
plot(v2a1_jskf_time,v2a1_jskf_ft,'--','LineWidth',2)
xlabel('time [s]')
ylabel('friction torque[Nm]')
xlim([0 inf])
title('Trajectory 4')
set(gca,'FontSize',20)

suptitle('Joint-3 friction torque')

%%

figure ('Name','estimated torque joint 3')

subplot(2,2,1)
hold on
plot(v1a1_org_time,v1a1_org_et)
plot(v1a1_jskf_time,v1a1_jskf_et,'--','LineWidth',2)
xlabel('time [s]')
ylabel('estimated torque[Nm]')
xlim([0 inf])
ylim([-inf inf])
title('Trajectory 1')
set(gca,'FontSize',20)

subplot(2,2,2)
hold on
plot(v2a2_org_time,v2a2_org_et)
plot(v2a2_jskf_time,v2a2_jskf_et,'--','LineWidth',2)
xlabel('time [s]')
ylabel('estimated torque[Nm]')
legend('raw signal','kalman filter')
xlim([0 inf])
ylim([-inf inf])
title('Trajectory 2')
set(gca,'FontSize',20)

subplot(2,2,3)
hold on
plot(v1a2_org_time,v1a2_org_et)
plot(v1a2_jskf_time,v1a2_jskf_et,'--','LineWidth',2)
xlabel('time [s]')
ylabel('estimated torque[Nm]')
xlim([0 inf])
ylim([-inf inf])
title('Trajectory 3')
set(gca,'FontSize',20)

subplot(2,2,4)
hold on
plot(v2a1_org_time,v2a1_org_et)
plot(v2a1_jskf_time,v2a1_jskf_et,'--','LineWidth',2)
xlabel('time [s]')
ylabel('estimated torque[Nm]')
xlim([0 inf])
ylim([-inf inf])
title('Trajectory 4')
set(gca,'FontSize',20)

suptitle('Joint-3 estimated torque')

%%
% 
% figure ('Name','friction torque joint 3')
% hold on
% plot(value{org,ft}.data(:,t),value{org,ft}.data(:,j3_t))
% plot(value{jskf,ft}.data(:,t),value{jskf,ft}.data(:,j3_t))
% xlabel('time [s]')
% ylabel('friction torque[Nm]')
% legend('raw signal','kalman filter')
% set(gca,'FontSize',16)
% suptitle('Joint-3 friction torque')
% 
% figure ('Name','estimated torque joint 3')
% hold on
% plot(value{org,et}.data(:,t),value{org,et}.data(:,j3_t))
% plot(value{jskf,et}.data(:,t),value{jskf,et}.data(:,j3_t))
% xlabel('time [s]')
% ylabel('estimated torque[Nm]')
% legend('raw signal','kalman filter')
% set(gca,'FontSize',16)
% suptitle('Joint-3 estimated torque')
