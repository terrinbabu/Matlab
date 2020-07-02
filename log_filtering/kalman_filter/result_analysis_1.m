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

path = '/home/terrin/projects/virtual_force_sensor/files/log_filtering/kalman_filter/0kg_all/';

value{org,it} = importdata(strcat(path,'ijt_0kg_all.txt'));
value{org,ft} = importdata(strcat(path,'fjt_0kg_all.txt'));
value{org,et} = importdata(strcat(path,'ejt_0kg_all.txt'));
value{org,mjs} = importdata(strcat(path,'mjt_0kg_all.txt'));
value{org,rt} = importdata(strcat(path,'rjt_0kg_all.txt'));
value{org,vt} = importdata(strcat(path,'vf_0kg_all.txt'));

value{jskf,it} = importdata(strcat(path,'ijt_0kg_all_jskf.txt'));
value{jskf,ft} = importdata(strcat(path,'fjt_0kg_all_jskf.txt'));
value{jskf,et} = importdata(strcat(path,'ejt_0kg_all_jskf.txt'));
value{jskf,mjs} = importdata(strcat(path,'mjt_0kg_all_jskf.txt'));
value{jskf,rt} = importdata(strcat(path,'rjt_0kg_all_jskf.txt'));
value{jskf,vt} = importdata(strcat(path,'vf_0kg_all_jskf.txt'));

value{pikf,it} = importdata(strcat(path,'ijt_0kg_all_pikf.txt'));
value{pikf,ft} = importdata(strcat(path,'fjt_0kg_all_pikf.txt'));
value{pikf,et} = importdata(strcat(path,'ejt_0kg_all_pikf.txt'));
value{pikf,mjs} = importdata(strcat(path,'mjt_0kg_all_pikf.txt'));
value{pikf,rt} = importdata(strcat(path,'rjt_0kg_all_pikf.txt'));
value{pikf,vt} = importdata(strcat(path,'vf_0kg_all_pikf.txt'));

value{pikf_jskf,it} = importdata(strcat(path,'ijt_0kg_all_pikf_jskf.txt'));
value{pikf_jskf,ft} = importdata(strcat(path,'fjt_0kg_all_pikf_jskf.txt'));
value{pikf_jskf,et} = importdata(strcat(path,'ejt_0kg_all_pikf_jskf.txt'));
value{pikf_jskf,mjs} = importdata(strcat(path,'mjt_0kg_all_pikf_jskf.txt'));
value{pikf_jskf,rt} = importdata(strcat(path,'rjt_0kg_all_pikf_jskf.txt'));
value{pikf_jskf,vt} = importdata(strcat(path,'vf_0kg_all_pikf_jskf.txt'));

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

%% position,Velocity and accerleration

figure ('Name','position joint 3')
hold on
plot(value{org,mjs}.data(:,t),value{org,mjs}.data(:,j3_q))
plot(value{jskf,mjs}.data(:,t),value{jskf,mjs}.data(:,j3_q))
legend('org','jskf')

figure ('Name','Velocity joint 3')
hold on
plot(value{org,mjs}.data(:,t),value{org,mjs}.data(:,j3_Dq))
plot(value{jskf,mjs}.data(:,t),value{jskf,mjs}.data(:,j3_Dq))
legend('org','jskf')

figure ('Name','accerleration joint 3')
hold on
plot(value{org,it}.data(:,t),value{org,it}.data(:,j3_DDq))
plot(value{jskf,it}.data(:,t),value{jskf,it}.data(:,j3_DDq))
legend('org','jskf')

figure ('Name','measured tau joint 3')
hold on
plot(value{org,mjs}.data(:,t),value{org,mjs}.data(:,j3_t))
plot(value{jskf,mjs}.data(:,t),value{jskf,mjs}.data(:,j3_t))
legend('org','jskf')

%% torques 

% inertial 

figure ('Name','inertial tau joint 3')
subplot(3,1,1)
hold on
plot(value{org,it}.data(:,t),value{org,it}.data(:,j3_t))
plot(value{jskf,it}.data(:,t),value{jskf,it}.data(:,j3_t))
legend('org','jskf')
subplot(3,1,2)
hold on
plot(value{org,it}.data(:,t),value{org,it}.data(:,j3_t))
plot(value{pikf,it}.data(:,t),value{pikf,it}.data(:,j3_t))
legend('org','pikf')
subplot(3,1,3)
hold on
plot(value{org,it}.data(:,t),value{org,it}.data(:,j3_t))
plot(value{pikf_jskf,it}.data(:,t),value{pikf_jskf,it}.data(:,j3_t))
legend('org','pikf+jskf')

% frictional

figure ('Name','friction tau joint 3')
subplot(3,1,1)
hold on
plot(value{org,ft}.data(:,t),value{org,ft}.data(:,j3_t))
plot(value{jskf,ft}.data(:,t),value{jskf,ft}.data(:,j3_t))
legend('org','jskf')
subplot(3,1,2)
hold on
plot(value{org,ft}.data(:,t),value{org,ft}.data(:,j3_t))
plot(value{pikf,ft}.data(:,t),value{pikf,ft}.data(:,j3_t))
legend('org','pikf')
subplot(3,1,3)
hold on
plot(value{org,ft}.data(:,t),value{org,ft}.data(:,j3_t))
plot(value{pikf_jskf,ft}.data(:,t),value{pikf_jskf,ft}.data(:,j3_t))
legend('org','pikf+jskf')

% estimated

figure ('Name','estimated tau joint 3')
subplot(3,1,1)
hold on
plot(value{org,et}.data(:,t),value{org,et}.data(:,j3_t))
plot(value{jskf,et}.data(:,t),value{jskf,et}.data(:,j3_t))
legend('org','jskf')
subplot(3,1,2)
hold on
plot(value{org,et}.data(:,t),value{org,et}.data(:,j3_t))
plot(value{pikf,et}.data(:,t),value{pikf,et}.data(:,j3_t))
legend('org','pikf')
subplot(3,1,3)
hold on
plot(value{org,et}.data(:,t),value{org,et}.data(:,j3_t))
plot(value{pikf_jskf,et}.data(:,t),value{pikf_jskf,et}.data(:,j3_t))
legend('pikf','pikf+jskf')

% residual 

figure ('Name','residual tau joint 3')
subplot(3,1,1)
hold on
plot(value{org,rt}.data(:,t),value{org,rt}.data(:,j3_t))
plot(value{jskf,rt}.data(:,t),value{jskf,rt}.data(:,j3_t))
legend('org','jskf')
subplot(3,1,2)
hold on
plot(value{org,rt}.data(:,t),value{org,rt}.data(:,j3_t))
plot(value{pikf,rt}.data(:,t),value{pikf,rt}.data(:,j3_t))
legend('org','pikf')
subplot(3,1,3)
hold on
plot(value{org,rt}.data(:,t),value{org,rt}.data(:,j3_t))
plot(value{pikf_jskf,rt}.data(:,t),value{pikf_jskf,rt}.data(:,j3_t))
legend('pikf','pikf+jskf')

% measured Vs estimated

figure ('Name','org measured Vs estimated tau - joint 3')
subplot(2,2,1)
hold on
plot(value{org,mjs}.data(:,t),value{org,mjs}.data(:,j3_t))
plot(value{org,et}.data(:,t),value{org,et}.data(:,j3_t))
legend('measured org','estimated org')
subplot(2,2,2)
hold on
plot(value{org,mjs}.data(:,t),value{org,mjs}.data(:,j3_t))
plot(value{jskf,et}.data(:,t),value{jskf,et}.data(:,j3_t))
legend('measured org','estimated jskf')
subplot(2,2,3)
hold on
plot(value{org,mjs}.data(:,t),value{org,mjs}.data(:,j3_t))
plot(value{pikf,et}.data(:,t),value{pikf,et}.data(:,j3_t))
legend('measured org','estimated pikf')
subplot(2,2,4)
hold on
plot(value{org,mjs}.data(:,t),value{org,mjs}.data(:,j3_t))
plot(value{pikf_jskf,et}.data(:,t),value{pikf_jskf,et}.data(:,j3_t))
legend('measured org','estimated pikf+jskf')

figure ('Name','jskf measured Vs estimated tau - joint 3')
subplot(2,2,1)
hold on
plot(value{jskf,mjs}.data(:,t),value{jskf,mjs}.data(:,j3_t))
plot(value{org,et}.data(:,t),value{org,et}.data(:,j3_t))
legend('measured jskf','estimated org')
subplot(2,2,2)
hold on
plot(value{jskf,mjs}.data(:,t),value{jskf,mjs}.data(:,j3_t))
plot(value{jskf,et}.data(:,t),value{jskf,et}.data(:,j3_t))
legend('measured jskf','estimated jskf')
subplot(2,2,3)
hold on
plot(value{jskf,mjs}.data(:,t),value{jskf,mjs}.data(:,j3_t))
plot(value{pikf,et}.data(:,t),value{pikf,et}.data(:,j3_t))
legend('measured jskf','estimated pikf')
subplot(2,2,4)
hold on
plot(value{jskf,mjs}.data(:,t),value{jskf,mjs}.data(:,j3_t))
plot(value{pikf_jskf,et}.data(:,t),value{pikf_jskf,et}.data(:,j3_t))
legend('measured jskf','estimated pikf+jskf')
