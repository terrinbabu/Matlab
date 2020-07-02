clc;
clear all;
close all;

%% time,position,velocity,torque,accelearation

seq = 1;
t = 2;

j1_q = 3; j1_Dq = 9;  j1_t = 15; j1_DDq = 21;
j2_q = 4; j2_Dq = 10; j2_t = 16; j2_DDq = 22;
j3_q = 5; j3_Dq = 11; j3_t = 17; j3_DDq = 23;
j4_q = 6; j4_Dq = 12; j4_t = 18; j4_DDq = 24;
j5_q = 7; j5_Dq = 13; j5_t = 19; j5_DDq = 25;
j6_q = 8; j6_Dq = 14; j6_t = 20; j6_DDq = 26;

%% input 

it_value = importdata('ijt_0kg_all.txt');

q = it_value.data(:,j3_q);
Dq = it_value.data(:,j3_Dq);
DDq  = it_value.data(:,j3_DDq);

%% write txt file

% fileID = fopen('Dq.txt','w');
% % fprintf(fileID,'%6s\n','Dq');
% fprintf(fileID,'%6.6f\n',Dq);
% fclose(fileID);

%% zero phase filter

% % % Moving average forward filter
% % 
% % N = 15;
% % MA_Dq = Dq;
% % size_Dq = size(Dq);
% % 
% % for n = N+1:size_Dq(1)
% %     
% %     MA_Dq(n)=0;
% %     
% %     for i = 0:N-1
% %         temp = Dq(n-i)/N;
% %         MA_Dq(n) = MA_Dq(n) + temp; 
% %     end
% % end 
% % 
% % % flip (flipud)
% % 
% % flip_MA_Dq = ones(size_Dq(1),1);
% % n=size_Dq(1);
% % 
% % for i = 1:size_Dq(1)
% %     flip_MA_Dq(n)=MA_Dq(i);
% %     n = n-1;
% % end
% % 
% % % Moving average backward filter
% % 
% % flip_FF_Dq = flip_MA_Dq;
% % 
% % 
% % 
% % for n = N+1:size_Dq(1)
% %     
% %     flip_FF_Dq(n)=0;
% %     
% %     for i = 0:N-1
% %         temp = flip_MA_Dq(n-i)/N;
% %         flip_FF_Dq(n) = flip_FF_Dq(n) + temp; 
% %     end
% % end
% % 
% % % flip (flipud)
% % 
% % FF_Dq = ones(size_Dq(1),1);
% % n=size_Dq(1);
% % 
% % for i = 1:size_Dq(1)
% %     FF_Dq(n)=flip_FF_Dq(i);
% %     n = n-1;
% % end
% % 
% % % results
% % 
% % figure('Name','filtered Dq');
% % hold on
% % plot(Dq);
% % plot(MA_Dq);
% % plot(FF_Dq);
% % xlabel('seq')
% % ylabel('Dq')
% % legend('Dq','ma-Dq','ff-Dq')

%% filtfilt filter

n=15;
coeffMA = ones(1,n)/n;

ff_q = filtfilt(coeffMA,1,q);
ff_Dq = filtfilt(coeffMA,1,Dq);
ff_DDq = filtfilt(coeffMA,1,DDq);

%% ploting the results

figure('Name','filtered Position');
hold on
plot(q);
plot(ff_q);
xlabel('seq')
ylabel('q')
legend('q','ff-q')

figure('Name','filtered Velocity');
hold on
plot(Dq);
plot(ff_Dq);
xlabel('seq')
ylabel('Dq')
legend('Dq','ff Dq')

figure('Name','filtered Accelearation');
hold on
plot(DDq);
plot(ff_DDq);
xlabel('seq')
ylabel('DDq')
legend('DDq','ff DDq')