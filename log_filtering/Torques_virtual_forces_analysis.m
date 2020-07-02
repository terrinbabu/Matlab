clc;
clear all;
close all;

%% input %%

% i - inertia, f- frictional, e - effective/estimated, m - motor/measured
% r- residual, vf - virtual force, fxxxx - filtered 

%%
roundn = @(x,n) round(x*10^n)./10^n;

it_value = importdata('ijt_no_load_sinjall_v5_1.txt');
ft_value = importdata('fjt_no_load_sinjall_v5_1.txt');
et_value = importdata('ejt_no_load_sinjall_v5_1.txt');
mjs_value = importdata('mjt_no_load_sinjall_v5_1.txt');
rt_value = importdata('rjt_no_load_sinjall_v5_1.txt');
vf_value = importdata('vf_no_load_sinjall_v5_1.txt');

start_trim = 0;     
end_trim = 0;

it_time = it_value.data(start_trim+1:end-end_trim,2);
it1 = roundn(it_value.data(start_trim+1:end-end_trim,15),3);
it2 = roundn(it_value.data(start_trim+1:end-end_trim,16),3);
it3 = roundn(it_value.data(start_trim+1:end-end_trim,17),3);
it4 = roundn(it_value.data(start_trim+1:end-end_trim,18),3);
it5 = roundn(it_value.data(start_trim+1:end-end_trim,19),3);
it6 = roundn(it_value.data(start_trim+1:end-end_trim,20),3);
iDDq1  = roundn(it_value.data(start_trim+1:end-end_trim,21),3);
iDDq2  = roundn(it_value.data(start_trim+1:end-end_trim,22),3);
iDDq3  = roundn(it_value.data(start_trim+1:end-end_trim,23),3);
iDDq4  = roundn(it_value.data(start_trim+1:end-end_trim,24),3);
iDDq5  = roundn(it_value.data(start_trim+1:end-end_trim,25),3);
iDDq6  = roundn(it_value.data(start_trim+1:end-end_trim,26),3);

ft_time = ft_value.data(start_trim+1:end-end_trim,2);
ft1 = roundn(ft_value.data(start_trim+1:end-end_trim,15),3);
ft2 = roundn(ft_value.data(start_trim+1:end-end_trim,16),3);
ft3 = roundn(ft_value.data(start_trim+1:end-end_trim,17),3);
ft4 = roundn(ft_value.data(start_trim+1:end-end_trim,18),3);
ft5 = roundn(ft_value.data(start_trim+1:end-end_trim,19),3);
ft6 = roundn(ft_value.data(start_trim+1:end-end_trim,20),3);

et_time = et_value.data(start_trim+1:end-end_trim,2);
et1 = roundn(et_value.data(start_trim+1:end-end_trim,15),3);
et2 = roundn(et_value.data(start_trim+1:end-end_trim,16),3);
et3 = roundn(et_value.data(start_trim+1:end-end_trim,17),3);
et4 = roundn(et_value.data(start_trim+1:end-end_trim,18),3);
et5 = roundn(et_value.data(start_trim+1:end-end_trim,19),3);
et6 = roundn(et_value.data(start_trim+1:end-end_trim,20),3);

mjs_seq = mjs_value.data(start_trim+1:end-end_trim,1);
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

rt_time = rt_value.data(start_trim+1:end-end_trim,2);
rt1 = roundn(rt_value.data(start_trim+1:end-end_trim,15),3);
rt2 = roundn(rt_value.data(start_trim+1:end-end_trim,16),3);
rt3 = roundn(rt_value.data(start_trim+1:end-end_trim,17),3);
rt4 = roundn(rt_value.data(start_trim+1:end-end_trim,18),3);
rt5 = roundn(rt_value.data(start_trim+1:end-end_trim,19),3);
rt6 = roundn(rt_value.data(start_trim+1:end-end_trim,20),3);

vf_time = vf_value.data(start_trim+1:end-end_trim,2);
fx = roundn(vf_value.data(start_trim+1:end-end_trim,3),3);
fy = roundn(vf_value.data(start_trim+1:end-end_trim,4),3);
fz = roundn(vf_value.data(start_trim+1:end-end_trim,5),3);

eval_it_time = it_time-it_time(1);
eval_ft_time = ft_time-ft_time(1);
eval_et_time = et_time-et_time(1);
eval_mjs_time = mjs_time-mjs_time(1);
eval_rt_time = rt_time-rt_time(1);
eval_vf_time = vf_time - vf_time(1);


%% Plots %%

%% Position 

figure ('Name','Position');

q_ylim_l = -4;
q_ylim_u = 4;

for i = 1:6
    subplot(3,2,i)
    plot(eval_mjs_time,roundn(mjs_value.data(start_trim+1:end-end_trim,i+2),3))
    xlabel('time')
    ylabel(strcat('q',int2str(i)))
    legend(strcat('q',int2str(i)))
%     ylim([q_ylim_l,q_ylim_u])
end

suptitle('Joint Postion')

%% Velocity

figure ('Name','Velocity')

Dq_ylim_l = -2;
Dq_ylim_u = 2;

for i = 1:6
    subplot(3,2,i)
    plot(eval_mjs_time,roundn(mjs_value.data(start_trim+1:end-end_trim,i+8),3))
    xlabel('time')
    ylabel(strcat('Dq',int2str(i)))
    legend(strcat('Dq',int2str(i)))
%     ylim([Dq_ylim_l,Dq_ylim_u])
end

suptitle('Joint Velocity')

%% Acceleration

figure ('Name','Acceleration')

DDq_ylim_l = -20;
DDq_ylim_u = 20;

for i = 1:6
    subplot(3,2,i)
    plot(eval_it_time,roundn(it_value.data(start_trim+1:end-end_trim,i+20),3))
    xlabel('time')
    ylabel(strcat('DDq',int2str(i)))
    legend(strcat('DDq',int2str(i)))
%     ylim([DDq_ylim_l,DDq_ylim_u])
end

suptitle('Joint Acceleration')


%% Torques

for i = 1:6

    figure ('Name',strcat('Torque - Joint ',int2str(i)))

    subplot(5,1,1)
    plot(eval_it_time,roundn(it_value.data(start_trim+1:end-end_trim,i+14),3))
    ylabel('inertial')

    subplot(5,1,2)
    plot(eval_ft_time,roundn(ft_value.data(start_trim+1:end-end_trim,i+14),3))
    ylabel('frictional')

    subplot(5,1,3)
    plot(eval_et_time,roundn(et_value.data(start_trim+1:end-end_trim,i+14),3))
    ylabel('estimated')

    subplot(5,1,4)
    plot(eval_mjs_time,roundn(mjs_value.data(start_trim+1:end-end_trim,i+14),3))
    ylabel('measured')

    subplot(5,1,5)
    plot(eval_rt_time,roundn(rt_value.data(start_trim+1:end-end_trim,i+14),3))
    ylabel('residual')

    for j = 1:5
        subplot(5,1,j)
        xlabel('time')
    end

    suptitle(strcat('Torque - Joint ',int2str(i)))

end


figure ('Name','Residual Torques of all Joints')

for i = 1:6
    subplot(3,2,i)
    plot(eval_rt_time,roundn(rt_value.data(start_trim+1:end-end_trim,i+14),3))
    xlabel('time')
    ylabel(strcat('rt',int2str(i)))
    legend(strcat('rt',int2str(i)))
end

suptitle('Residual Torques of all Joints')

%% Virtual Force

vf_ylim_l = -400;
vf_ylim_u = 400;

figure ('Name','Virtual Force')

subplot(3,1,1)
plot(eval_vf_time,fx)
ylabel('fx')

subplot(3,1,2)
plot(eval_vf_time,fy)
ylabel('fy')

subplot(3,1,3)
plot(eval_vf_time,fz)
ylabel('fz')

for i = 1:3
    subplot(3,1,i)
    xlabel('time')
%     ylim([vf_ylim_l,vf_ylim_u])
end

suptitle('Virtual Force')

%% Filtered joint states

%% Input 

fit_value = importdata('ijt_filtered_no_load_sinjall_v5_1_DDqfromDq.txt');
fft_value = importdata('fjt_filtered_no_load_sinjall_v5_1_DDqfromDq.txt');
fet_value = importdata('ejt_filtered_no_load_sinjall_v5_1_DDqfromDq.txt');
fmjs_value = importdata('mjt_filtered_no_load_sinjall_v5_1_DDqfromDq.txt');
frt_value = importdata('rjt_filtered_no_load_sinjall_v5_1_DDqfromDq.txt');
fvf_value = importdata('vf_filtered_no_load_sinjall_v5_1_DDqfromDq.txt');

start_trim = 0; 
end_trim = 0;

fit_time = fit_value.data(start_trim+1:end-end_trim,2);
fit1 = roundn(fit_value.data(start_trim+1:end-end_trim,15),3);
fit2 = roundn(fit_value.data(start_trim+1:end-end_trim,16),3);
fit3 = roundn(fit_value.data(start_trim+1:end-end_trim,17),3);
fit4 = roundn(fit_value.data(start_trim+1:end-end_trim,18),3);
fit5 = roundn(fit_value.data(start_trim+1:end-end_trim,19),3);
fit6 = roundn(fit_value.data(start_trim+1:end-end_trim,20),3);
fiDDq1  = roundn(fit_value.data(start_trim+1:end-end_trim,21),3);
fiDDq2  = roundn(fit_value.data(start_trim+1:end-end_trim,22),3);
fiDDq3  = roundn(fit_value.data(start_trim+1:end-end_trim,23),3);
fiDDq4  = roundn(fit_value.data(start_trim+1:end-end_trim,24),3);
fiDDq5  = roundn(fit_value.data(start_trim+1:end-end_trim,25),3);
fiDDq6  = roundn(fit_value.data(start_trim+1:end-end_trim,26),3);

fft_time = fft_value.data(start_trim+1:end-end_trim,2);
fft1 = roundn(fft_value.data(start_trim+1:end-end_trim,15),3);
fft2 = roundn(fft_value.data(start_trim+1:end-end_trim,16),3);
fft3 = roundn(fft_value.data(start_trim+1:end-end_trim,17),3);
fft4 = roundn(fft_value.data(start_trim+1:end-end_trim,18),3);
fft5 = roundn(fft_value.data(start_trim+1:end-end_trim,19),3);
fft6 = roundn(fft_value.data(start_trim+1:end-end_trim,20),3);

fet_time = fet_value.data(start_trim+1:end-end_trim,2);
fet1 = roundn(fet_value.data(start_trim+1:end-end_trim,15),3);
fet2 = roundn(fet_value.data(start_trim+1:end-end_trim,16),3);
fet3 = roundn(fet_value.data(start_trim+1:end-end_trim,17),3);
fet4 = roundn(fet_value.data(start_trim+1:end-end_trim,18),3);
fet5 = roundn(fet_value.data(start_trim+1:end-end_trim,19),3);
fet6 = roundn(fet_value.data(start_trim+1:end-end_trim,20),3);

fmjs_seq = fmjs_value.data(start_trim+1:end-end_trim,1);
fmjs_time = fmjs_value.data(start_trim+1:end-end_trim,2);
fmq1 = roundn(fmjs_value.data(start_trim+1:end-end_trim,3),3);
fmq2 = roundn(fmjs_value.data(start_trim+1:end-end_trim,4),3);
fmq3 = roundn(fmjs_value.data(start_trim+1:end-end_trim,5),3);
fmq4 = roundn(fmjs_value.data(start_trim+1:end-end_trim,6),3);
fmq5 = roundn(fmjs_value.data(start_trim+1:end-end_trim,7),3);
fmq6 = roundn(fmjs_value.data(start_trim+1:end-end_trim,8),3);
fmDq1 = roundn(fmjs_value.data(start_trim+1:end-end_trim,9),3);
fmDq2 = roundn(fmjs_value.data(start_trim+1:end-end_trim,10),3);
fmDq3 = roundn(fmjs_value.data(start_trim+1:end-end_trim,11),3);
fmDq4 = roundn(fmjs_value.data(start_trim+1:end-end_trim,12),3);
fmDq5 = roundn(fmjs_value.data(start_trim+1:end-end_trim,13),3);
fmDq6 = roundn(fmjs_value.data(start_trim+1:end-end_trim,14),3);
fmt1 = roundn(fmjs_value.data(start_trim+1:end-end_trim,15),3);
fmt2 = roundn(fmjs_value.data(start_trim+1:end-end_trim,16),3);
fmt3 = roundn(fmjs_value.data(start_trim+1:end-end_trim,17),3);
fmt4 = roundn(fmjs_value.data(start_trim+1:end-end_trim,18),3);
fmt5 = roundn(fmjs_value.data(start_trim+1:end-end_trim,19),3);
fmt6 = roundn(fmjs_value.data(start_trim+1:end-end_trim,20),3);

frt_time = frt_value.data(start_trim+1:end-end_trim,2);
frt1 = roundn(frt_value.data(start_trim+1:end-end_trim,15),3);
frt2 = roundn(frt_value.data(start_trim+1:end-end_trim,16),3);
frt3 = roundn(frt_value.data(start_trim+1:end-end_trim,17),3);
frt4 = roundn(frt_value.data(start_trim+1:end-end_trim,18),3);
frt5 = roundn(frt_value.data(start_trim+1:end-end_trim,19),3);
frt6 = roundn(frt_value.data(start_trim+1:end-end_trim,20),3);

fvf_time = fvf_value.data(start_trim+1:end-end_trim,2);
ffx = roundn(fvf_value.data(start_trim+1:end-end_trim,3),3);
ffy = roundn(fvf_value.data(start_trim+1:end-end_trim,4),3);
ffz = roundn(fvf_value.data(start_trim+1:end-end_trim,5),3);

eval_fit_time = fit_time - fit_time(1);
eval_fft_time = fft_time - fft_time(1);
eval_fet_time = fet_time - fet_time(1);
eval_fmjs_time = fmjs_time - fmjs_time(1);
eval_frt_time = frt_time - frt_time(1);
eval_fvf_time = fvf_time - fvf_time(1);

%% Position 

figure ('Name','Position with Filter');

for i = 1:6
    
    subplot(3,2,i)
    plot(eval_mjs_time,roundn(mjs_value.data(start_trim+1:end-end_trim,i+2),3))
    hold on
    plot(eval_fmjs_time,roundn(fmjs_value.data(start_trim+1:end-end_trim,i+2),3))    
    
    xlabel('time')
    ylabel(strcat('q',int2str(i),'/fq',int2str(i)))
    legend(strcat('q',int2str(i)),strcat('fq',int2str(i)))
end

suptitle('Joint Postion with Filter')

%% Velocity

figure ('Name','Velocity  with Filter')

for i = 1:6
    
    subplot(3,2,i)
    plot(eval_mjs_time,roundn(mjs_value.data(start_trim+1:end-end_trim,i+8),3))
    hold on
    plot(eval_fmjs_time,roundn(fmjs_value.data(start_trim+1:end-end_trim,i+8),3))
    
    xlabel('time')
    ylabel(strcat('Dq',int2str(i),'/fDq',int2str(i)))
    legend(strcat('Dq',int2str(i)),strcat('fDq',int2str(i)))
end

suptitle('Joint Velocity  with Filter')

%% Acceleration

figure ('Name','Acceleration  with Filter')

for i = 1:6
    
    subplot(3,2,i)
    plot(eval_it_time,roundn(it_value.data(start_trim+1:end-end_trim,i+20),3))
    hold on
    plot(eval_fit_time,roundn(fit_value.data(start_trim+1:end-end_trim,i+20),3))
    
    xlabel('time')
    ylabel(strcat('DDq',int2str(i),'/fDDq',int2str(i)))
    legend(strcat('DDq',int2str(i)),strcat('fDDq',int2str(i)))
end

suptitle('Joint Acceleration  with Filter')

%% 

figure ('Name','Inertial torque with filter')

for i = 1:6
    
    subplot(3,2,i)
    plot(eval_it_time,roundn(it_value.data(start_trim+1:end-end_trim,i+14),3))
    hold on
    plot(eval_fit_time,roundn(fit_value.data(start_trim+1:end-end_trim,i+14),3))
    
    xlabel('time')
    ylabel(strcat('it',int2str(i),'/fit',int2str(i)))
    legend(strcat('it',int2str(i)),strcat('fit',int2str(i)))
end

suptitle('Inertial torque with filter')

%%
figure ('Name','Frictional torque with filter')

for i = 1:6
    
    subplot(3,2,i)
    plot(eval_ft_time,roundn(ft_value.data(start_trim+1:end-end_trim,i+14),3))
    hold on
    plot(eval_fft_time,roundn(fft_value.data(start_trim+1:end-end_trim,i+14),3))
    
    xlabel('time')
    ylabel(strcat('ft',int2str(i),'/fft',int2str(i)))
    legend(strcat('ft',int2str(i)),strcat('fft',int2str(i)))
end

suptitle('Frictional torque with filter')

%%
figure ('Name','Estimated torque with filter')

for i = 1:6
    
    subplot(3,2,i)
    plot(eval_et_time,roundn(et_value.data(start_trim+1:end-end_trim,i+14),3))
    hold on
    plot(eval_fet_time,roundn(fet_value.data(start_trim+1:end-end_trim,i+14),3))
    
    xlabel('time')
    ylabel(strcat('et',int2str(i),'/fet',int2str(i)))
    legend(strcat('et',int2str(i)),strcat('fet',int2str(i)))
end

suptitle('Estimated torque with filter')

%%

figure ('Name','Measured torque with filter')

for i = 1:6
    
    subplot(3,2,i)
    plot(eval_mjs_time,roundn(mjs_value.data(start_trim+1:end-end_trim,i+14),3))
    hold on
    plot(eval_fmjs_time,roundn(fmjs_value.data(start_trim+1:end-end_trim,i+14),3))
    
    xlabel('time')
    ylabel(strcat('mt',int2str(i),'/fmt',int2str(i)))
    legend(strcat('mt',int2str(i)),strcat('fmt',int2str(i)))
end

suptitle('Measured torque with filter')

%%

figure ('Name','Residual torque with filter')

for i = 1:6
    
    subplot(3,2,i)
    plot(eval_rt_time,roundn(rt_value.data(start_trim+1:end-end_trim,i+14),3))
    hold on
    plot(eval_frt_time,roundn(frt_value.data(start_trim+1:end-end_trim,i+14),3))
    
    xlabel('time')
    ylabel(strcat('rt',int2str(i),'/frt',int2str(i)))
    legend(strcat('rt',int2str(i)),strcat('frt',int2str(i)))
end

suptitle('Residual torque with filter')

%%

figure ('Name','Virtual Force with filter')

subplot(3,1,1)
plot(eval_vf_time,fx)
hold on
plot(eval_fvf_time,ffx)
ylabel('fx')
legend('fx','ffx')

subplot(3,1,2)
plot(eval_vf_time,fy)
hold on
plot(eval_fvf_time,ffy)
ylabel('fy')
legend('fy','ffy')

subplot(3,1,3)
plot(eval_vf_time,fz)
hold on
plot(eval_fvf_time,ffz)
ylabel('fz')
legend('fz','ffz')

for i = 1:3
    subplot(3,1,i)
    xlabel('time')
end

suptitle('Virtual Force with filter')

%% Figure for Presentation 

% % % %% Position
% % % 
% % % figure ('Name','Position');
% % % 
% % % q_ylim_l = -4;
% % % q_ylim_u = 4;
% % % 
% % % for i = 1:6
% % %     subplot(3,2,i)
% % %     plot(eval_mjs_time,roundn(mjs_value.data(start_trim+1:end-end_trim,i+2),3))
% % %     xlabel('time')
% % %     ylabel({strcat('Joint-',int2str(i));'Position [rad]'})
% % % %     legend(strcat('q',int2str(i)))
% % %     set(gca,'FontSize',16)
% % % end
% % % 
% % % suptitle('Joint Postion')
% % % 
% % % %% Velocity
% % % 
% % % figure ('Name','Velocity')
% % % 
% % % Dq_ylim_l = -2;
% % % Dq_ylim_u = 2;
% % % 
% % % for i = 1:6
% % %     subplot(3,2,i)
% % %     plot(eval_mjs_time,roundn(mjs_value.data(start_trim+1:end-end_trim,i+8),3))
% % %     xlabel('time')
% % %     ylabel({strcat('Joint-',int2str(i));'Velocity [rad/s]'})
% % % %     legend(strcat('Dq',int2str(i)))
% % %     set(gca,'FontSize',16)
% % % end
% % % 
% % % suptitle('Joint Velocity')
% % % 
% % % %% Acceleration
% % % 
% % % figure ('Name','Acceleration')
% % % 
% % % DDq_ylim_l = -20;
% % % DDq_ylim_u = 20;
% % % 
% % % for i = 1:6
% % %     subplot(3,2,i)
% % %     plot(eval_it_time,roundn(it_value.data(start_trim+1:end-end_trim,i+20),3))
% % %     xlabel('time')
% % %     ylabel({strcat('Joint-',int2str(i));'Acceleration [rad/s^2]'})
% % % %     legend(strcat('DDq',int2str(i)))
% % %     set(gca,'FontSize',16)
% % % %     ylim([DDq_ylim_l,DDq_ylim_u])
% % % end
% % % 
% % % suptitle('Joint Acceleration')
% % % set(gca,'FontSize',16)
% % % 
% % % 
% % %  
% % % figure ('Name','Virtual Force with filter')
% % % 
% % % start_trim = 500;
% % % end_trim = length(fx)- 1500;
% % % 
% % % fx = roundn(virtual_force_value.data(start_trim+1:end-end_trim,3),3);
% % % fy = roundn(virtual_force_value.data(start_trim+1:end-end_trim,4),3);
% % % fz = roundn(virtual_force_value.data(start_trim+1:end-end_trim,5),3);
% % % ffx = roundn(fvirtual_force_value.data(start_trim+1:end-end_trim,3),3);
% % % ffy = roundn(fvirtual_force_value.data(start_trim+1:end-end_trim,4),3);
% % % ffz = roundn(fvirtual_force_value.data(start_trim+1:end-end_trim,5),3);
% % % 
% % % 
% % % vf_ylim_l = -400;
% % % vf_ylim_u = 400;
% % % 
% % % subplot(3,1,1)
% % % plot(fx)
% % % hold on
% % % plot(ffx)
% % % ylabel('force x-axis (N)')
% % % legend('fx','joint filtered fx')
% % % set(gca,'FontSize',16)
% % % 
% % % subplot(3,1,2)
% % % plot(fy)
% % % hold on
% % % plot(ffy)
% % % ylabel('force y-axis (N)')
% % % legend('fy','joint filtered fy')
% % % set(gca,'FontSize',16)
% % % 
% % % subplot(3,1,3)
% % % plot(fz)
% % % hold on
% % % plot(ffz)
% % % ylabel('force z-axis (N)')
% % % legend('fz','joint filtered fz')
% % % set(gca,'FontSize',16)
% % % 
% % % for i = 1:3
% % %     subplot(3,1,i)
% % %     xlabel('time')
% % %     ylim([vf_ylim_l,vf_ylim_u])
% % % end
% % % 
% % % 
% % % 
% % % suptitle('Virtual Force')
% % % set(gca,'FontSize',16)
% % % 
% % % % % 
% % % % % 
% % % % % figure ('Name','Joint-2 Torques')
% % % % % 
% % % % % q_ylim_l = -2;
% % % % % q_ylim_u = 0.5;
% % % % % 
% % % % % 
% % % % % subplot(4,1,1)
% % % % % yyaxis left
% % % % % plot(eval_ijs_time,mq2);
% % % % % xlabel('time')
% % % % % ylabel('Position [rad]')
% % % % % ylim([q_ylim_l,q_ylim_u])
% % % % % 
% % % % % yyaxis right
% % % % % plot(eval_ijs_time,mDq2);
% % % % % set(gca,'FontSize',14)
% % % % % ylabel('Velocity [rad/s]')
% % % % % legend('Position ','Velocity')
% % % % % 
% % % % % 
% % % % % subplot(4,1,2)
% % % % % yyaxis left
% % % % % plot(eval_ijs_time,ieff2);
% % % % % xlabel('time')
% % % % % ylabel({'Inertial';'Torque [Nm]'})
% % % % % 
% % % % % yyaxis right
% % % % % plot(eval_ijs_time,feff2);
% % % % % set(gca,'FontSize',14)
% % % % % ylabel({'Friction';'Torque [Nm]'})
% % % % % legend('Inertial Torque','Friction Torque')
% % % % % 
% % % % % subplot(4,1,3)
% % % % % plot(eval_ijs_time,[eeff2,meff2]);
% % % % % xlabel('time')
% % % % % legend('Estimated Torque','Motor Torque')
% % % % % ylabel({'Estimated/Motor';'Torque [Nm]'})
% % % % % set(gca,'FontSize',14)
% % % % % 
% % % % % subplot(4,1,4)
% % % % % plot(eval_ijs_time,reff2)
% % % % % ylabel({'Residual';'Torque [Nm]'})
% % % % % set(gca,'FontSize',14)
% % % % % 
% % % % % for i = 1:4
% % % % %     subplot(4,1,i)
% % % % %     xlabel('time')
% % % % % end
% % % % % 
% % % % % % h = suptitle('Joint-2 Residual Torques Estimation');
% % % % % % set(h,'FontSize',20)
