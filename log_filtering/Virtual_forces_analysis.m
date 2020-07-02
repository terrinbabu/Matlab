% input

close all; clc;
roundn = @(x,n) round(x*10^n)./10^n;

joint_states_value = importdata('js_v10_-y.txt');
fjoint_states_value = importdata('ejs_v10_-y.txt');
%fjoint_states_value.data = joint_states_value.data;
virtual_force_value = importdata('vf_v10_-y.txt');
%sensor_force_value = importdata('sf_-z_kf_vel50_P1.txt');

start_trim = 0;
end_trim = 0;

js_seq = joint_states_value.data(start_trim+1:end-end_trim,1);
js_time = joint_states_value.data(start_trim+1:end-end_trim,2);
q0 = roundn(joint_states_value.data(start_trim+1:end-end_trim,3),3);
q1 = roundn(joint_states_value.data(start_trim+1:end-end_trim,4),3);
q2 = roundn(joint_states_value.data(start_trim+1:end-end_trim,5),3);
q3 = roundn(joint_states_value.data(start_trim+1:end-end_trim,6),3);
q4 = roundn(joint_states_value.data(start_trim+1:end-end_trim,7),3);
q5 = roundn(joint_states_value.data(start_trim+1:end-end_trim,8),3);
Dq0 = roundn(joint_states_value.data(start_trim+1:end-end_trim,9),3);
Dq1 = roundn(joint_states_value.data(start_trim+1:end-end_trim,10),3);
Dq2 = roundn(joint_states_value.data(start_trim+1:end-end_trim,11),3);
Dq3 = roundn(joint_states_value.data(start_trim+1:end-end_trim,12),3);
Dq4 = roundn(joint_states_value.data(start_trim+1:end-end_trim,13),3);
Dq5 = roundn(joint_states_value.data(start_trim+1:end-end_trim,14),3);
eff0 = roundn(joint_states_value.data(start_trim+1:end-end_trim,15),3);
eff1 = roundn(joint_states_value.data(start_trim+1:end-end_trim,16),3);
eff2 = roundn(joint_states_value.data(start_trim+1:end-end_trim,17),3);
eff3 = roundn(joint_states_value.data(start_trim+1:end-end_trim,18),3);
eff4 = roundn(joint_states_value.data(start_trim+1:end-end_trim,19),3);
eff5 = roundn(joint_states_value.data(start_trim+1:end-end_trim,20),3);

fjs_seq = fjoint_states_value.data(start_trim+1:end-end_trim,1);
fjs_time = fjoint_states_value.data(start_trim+1:end-end_trim,2);
fq0 = roundn(fjoint_states_value.data(start_trim+1:end-end_trim,3),3);
fq1 = roundn(fjoint_states_value.data(start_trim+1:end-end_trim,4),3);
fq2 = roundn(fjoint_states_value.data(start_trim+1:end-end_trim,5),3);
fq3 = roundn(fjoint_states_value.data(start_trim+1:end-end_trim,6),3);
fq4 = roundn(fjoint_states_value.data(start_trim+1:end-end_trim,7),3);
fq5 = roundn(fjoint_states_value.data(start_trim+1:end-end_trim,8),3);
fDq0 = roundn(fjoint_states_value.data(start_trim+1:end-end_trim,9),3);
fDq1 = roundn(fjoint_states_value.data(start_trim+1:end-end_trim,10),3);
fDq2 = roundn(fjoint_states_value.data(start_trim+1:end-end_trim,11),3);
fDq3 = roundn(fjoint_states_value.data(start_trim+1:end-end_trim,12),3);
fDq4 = roundn(fjoint_states_value.data(start_trim+1:end-end_trim,13),3);
fDq5 = roundn(fjoint_states_value.data(start_trim+1:end-end_trim,14),3);
feff0 = roundn(fjoint_states_value.data(start_trim+1:end-end_trim,15),3);
feff1 = roundn(fjoint_states_value.data(start_trim+1:end-end_trim,16),3);
feff2 = roundn(fjoint_states_value.data(start_trim+1:end-end_trim,17),3);
feff3 = roundn(fjoint_states_value.data(start_trim+1:end-end_trim,18),3);
feff4 = roundn(fjoint_states_value.data(start_trim+1:end-end_trim,19),3);
feff5 = roundn(fjoint_states_value.data(start_trim+1:end-end_trim,20),3);

vf_seq = virtual_force_value.data(start_trim+1:end-end_trim,1);
vf_time = virtual_force_value.data(start_trim+1:end-end_trim,2);
fx = roundn(virtual_force_value.data(start_trim+1:end-end_trim,3),3);
fy = roundn(virtual_force_value.data(start_trim+1:end-end_trim,4),3);
fz = roundn(virtual_force_value.data(start_trim+1:end-end_trim,5),3);
tx = roundn(virtual_force_value.data(start_trim+1:end-end_trim,6),3);
ty = roundn(virtual_force_value.data(start_trim+1:end-end_trim,7),3);
tz = roundn(virtual_force_value.data(start_trim+1:end-end_trim,8),3);

% sf_seq = sensor_force_value.data(start_trim+1:end-end_trim,1);
% sf_time = sensor_force_value.data(start_trim+1:end-end_trim,2);
% sfx = roundn(sensor_force_value.data(start_trim+1:end-end_trim,3),3);
% sfy = roundn(sensor_force_value.data(start_trim+1:end-end_trim,4),3);
% sfz = roundn(sensor_force_value.data(start_trim+1:end-end_trim,5),3);
% stx = roundn(sensor_force_value.data(start_trim+1:end-end_trim,6),3);
% sty = roundn(sensor_force_value.data(start_trim+1:end-end_trim,7),3);
% stz = roundn(sensor_force_value.data(start_trim+1:end-end_trim,8),3);

eval_js_time = js_time-js_time(1);
eval_fjs_time = fjs_time-fjs_time(1);
eval_vf_time = vf_time-vf_time(1);
% eval_sf_time = sf_time-sf_time(1);

%% frequency 

% % start_freq_check = 1;
% % 
% % for i = start_freq_check:(length(js_time)-1)
% %     freq_js(i) = 1/(js_time(i+1)-js_time(i));
% % end
% % 
% % for i = start_freq_check:(length(fjs_time)-1)
% %     freq_fjs(i) = 1/(fjs_time(i+1)-fjs_time(i));
% % end
% % 
% % for i = start_freq_check:(length(vf_time)-1)
% %     freq_vf(i) = 1/(vf_time(i+1)-vf_time(i));
% % end
% % 
% % for i = start_freq_check:(length(sf_time)-1)
% %     freq_sf(i) = 1/(sf_time(i+1)-sf_time(i));
% % end
% % 
% %  nSamp = 10;
% % MA_freq_js = MA_filter(freq_js,nSamp);
% % MA_freq_fjs = MA_filter(freq_fjs,nSamp);
% % MA_freq_vf = MA_filter(freq_vf,nSamp);
% % MA_freq_sf = MA_filter(freq_sf,nSamp);
% % 
% % plot_hz = figure ('Name','freq');
% % 
% % subplot(2,2,1)
% % plot(MA_freq_js)
% % ylabel('freq js')
% % 
% % for i = 1:4
% %     subplot(2,2,i)
% %     xlabel('time')
% % end
% % 
% % 
% % subplot(2,2,2)
% % plot(MA_freq_fjs)
% % ylabel('freq fjs')
% % 
% % subplot(2,2,3)
% % plot(MA_freq_vf)
% % ylabel('freq vf')
% % 
% % subplot(2,2,4)
% % plot(MA_freq_sf)
% % ylabel('freq sf')
% % 
% % for i = 1:4
% %     subplot(2,2,i)
% %     xlabel('time')
% % end


%% filtering %%

% req_signal = eff0;
% req_signal_time = js_time;
% freq_signal = feff0;
% freq_signal_time = fjs_time;
% 
% % MA Filter
% nSamp = 10;
% MA_filtered_signal = MA_filter(req_signal,nSamp);
% 
% FIR Filter
% FIR_order = 32;
% cut_freq = 20;
% FIR_filtered_signal = LP_FIR_filter(req_signal,req_signal_time,cut_freq,FIR_order);
% 
% % IIR Filter
% IIR_order = 2;
% cut_freq = 20;
% IIR_filtered_signal = LP_butter_filter(req_signal,req_signal_time,cut_freq,IIR_order);

%% Plots - filtering %%

% Power Spectral Density
%
% PSD(req_signal,req_signal_time);
% PSD(freq_signal,freq_signal_time);
% PSD(MA_filtered_signal,req_signal_time);
% PSD(FIR_filtered_signal,req_signal_time);
% PSD(IIR_filtered_signal,req_signal_time);
% 
% figure ('Name',getVarName(Dq0))
% 
% subplot(2,2,1)
% plot(req_signal_time,req_signal,freq_signal_time,freq_signal)
% xlabel('time') 
% ylabel('signal/filtered signal')
% legend('req signal','pre filtered signal')
% 
% subplot(2,2,2)
% plot(req_signal_time,[req_signal,MA_filtered_signal])
% xlabel('time') 
% ylabel('signal/filtered signal')
% legend('req signal','MA filtered signal')
% 
% subplot(2,2,3)
% plot(req_signal_time,[req_signal,FIR_filtered_signal])
% xlabel('time') 
% ylabel('signal/filtered signal')
% legend('req signal','FIR filtered signal')
% 
% subplot(2,2,4)
% plot(req_signal_time,[req_signal,IIR_filtered_signal])
% xlabel('time') 
% ylabel('signal/filtered signal')
% legend('req signal','IIR filtered signal')
% suptitle(getVarName(Dq0))

%% Plots - Joint state Analysis %%
 
q_ylim_l = -4;
q_ylim_u = 4;
Dq_ylim_l = -10;
Dq_ylim_u = 10;
eff_ylim_l = -150;
eff_ylim_u = 150;

plot_q = figure ('Name','Position');

subplot(3,2,1)
plot(eval_js_time,q0,'--')
hold on
plot(eval_fjs_time,fq0,':')

subplot(3,2,2)
plot(eval_js_time,q1,'--')
hold on
plot(eval_fjs_time,fq1,':')

subplot(3,2,3)
plot(eval_js_time,q2,'--')
hold on
plot(eval_fjs_time,fq2,':')

subplot(3,2,4)
plot(eval_js_time,q3,'--')
hold on
plot(eval_fjs_time,fq3,':')

subplot(3,2,5)
plot(eval_js_time,q4,'--')
hold on
plot(eval_fjs_time,fq4,':')

subplot(3,2,6)
plot(eval_js_time,q5,'--')
hold on
plot(eval_fjs_time,fq5,':')

for i = 1:6
    subplot(3,2,i)
    xlabel('time')
    ylabel(strcat('q',int2str(i-1),'/fq',int2str(i-1)))
    legend(strcat('q',int2str(i-1)),strcat('fq',int2str(i-1)))
    ylim([q_ylim_l,q_ylim_u])
end

suptitle('Joint Postion before and after filtering')

figure ('Name','Velocity')

subplot(3,2,1)
plot(eval_js_time,Dq0,'--')
hold on
plot(eval_fjs_time,fDq0,':')

subplot(3,2,2)
plot(eval_js_time,Dq1,'--')
hold on
plot(eval_fjs_time,fDq1,':')

subplot(3,2,3)
plot(eval_js_time,Dq2,'--')
hold on
plot(eval_fjs_time,fDq2,':')

subplot(3,2,4)
plot(eval_js_time,Dq3,'--')
hold on
plot(eval_fjs_time,fDq3,':')

subplot(3,2,5)
plot(eval_js_time,Dq4,'--')
hold on
plot(eval_fjs_time,fDq4,':')

subplot(3,2,6)
plot(eval_js_time,Dq5,'--')
hold on
plot(eval_fjs_time,fDq5,':')

for i = 1:6
    subplot(3,2,i)
    xlabel('time')
    ylabel(strcat('Dq',int2str(i-1),'/fDq',int2str(i-1)))
    legend(strcat('Dq',int2str(i-1)),strcat('fDq',int2str(i-1)))
    %ylim([Dq_ylim_l,Dq_ylim_u])
end

suptitle('Joint Velocity before and after filtering')

figure ('Name','Effort')

subplot(3,2,1)
plot(eval_js_time,eff0,'--')
hold on
plot(eval_fjs_time,feff0,':')

subplot(3,2,2)
plot(eval_js_time,eff1,'--')
hold on
plot(eval_fjs_time,feff1,':')

subplot(3,2,3)
plot(eval_js_time,eff2,'--')
hold on
plot(eval_fjs_time,feff2,':')

subplot(3,2,4)
plot(eval_js_time,eff3,'--')
hold on
plot(eval_fjs_time,feff3,':')

subplot(3,2,5)
plot(eval_js_time,eff4,'--')
hold on
plot(eval_fjs_time,feff4,':')

subplot(3,2,6)
plot(eval_js_time,eff5,'--')
hold on
plot(eval_fjs_time,feff5,':')

for i = 1:6
    subplot(3,2,i)
    xlabel('time')
    ylabel(strcat('eff',int2str(i-1),'/feff',int2str(i-1)))
    legend(strcat('eff',int2str(i-1)),strcat('feff',int2str(i-1)))
    %ylim([eff_ylim_l,eff_ylim_u])
end

suptitle('Joint Effort before and after filtering')

%% force analysis 

 vf_ylim_l = -300;
 vf_ylim_u = 300;


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
    %ylim([vf_ylim_l,vf_ylim_u])
end

suptitle('Virtual Force')

% % figure ('Name','Sensor Force')
% % 
% % subplot(3,1,1)
% % plot(eval_sf_time,sfx)
% % ylabel('sfx')
% % 
% % subplot(3,1,2)
% % plot(eval_sf_time,sfy)
% % ylabel('sfy')
% % 
% % subplot(3,1,3)
% % plot(eval_sf_time,sfz)
% % ylabel('sfz')
% % 
% % for i = 1:3
% %     subplot(3,1,i)
% %     xlabel('time')
% % %    ylim([vf_ylim_l,vf_ylim_u])
% % end
% % 
% % suptitle('Sensor Force')
% % 
% % figure ('Name','Sensor Force Vs Virtual Force')
% % 
% % subplot(3,1,1)
% % plot(eval_vf_time,fz*-1,'b')
% % hold on
% % plot(eval_sf_time,sfx,'g')
% % ylabel('sfx/vfz')
% % legend('-vfz','sfx')
% % 
% % subplot(3,1,2)
% % plot(eval_vf_time,fy,'b')
% % hold on
% % plot(eval_sf_time,sfy,'g')
% % ylabel('sfy/vfy')
% % legend('vfy','sfy')
% % 
% % subplot(3,1,3)
% % plot(eval_vf_time,fx,'b')
% % hold on
% % plot(eval_sf_time,sfz,'g')
% % ylabel('sfz/vfx')
% % legend('vfx','sfz')
% % 
% % for i = 1:3
% %     subplot(3,1,i)
% %     xlabel('time')
% % %    ylim([vf_ylim_l,vf_ylim_u])
% % end
% % 
% % suptitle('Sensor Force Vs Virtual Force')

%% creating subset
% % % 
% % % %100% override
subset_start_time = 6;
subset_end_time = 16;

vf_time_trim_1 = eval_vf_time(eval_vf_time>subset_start_time) ;
vf_start_trim_value = vf_time_trim_1(1);
vf_start_trim_index = find( eval_vf_time== vf_start_trim_value);

vf_time_trim_2 = eval_vf_time(eval_vf_time<subset_end_time);
vf_end_trim_value = vf_time_trim_2(end);
vf_end_trim_index = find( eval_vf_time== vf_end_trim_value);
% % % 
% % % sf_time_trim_1 = eval_sf_time(eval_sf_time>subset_start_time) ;
% % % sf_start_trim_value = sf_time_trim_1(1);
% % % sf_start_trim_index = find( eval_sf_time == sf_start_trim_value);
% % % 
% % % sf_time_trim_2 = eval_sf_time(eval_sf_time<subset_end_time);
% % % sf_end_trim_value = sf_time_trim_2(end);
% % % sf_end_trim_index = find( eval_sf_time == sf_end_trim_value);
% % % 
vf_start_trim_index = vf_start_trim_index(1);
vf_end_trim_index = vf_end_trim_index(1);
% % % sf_start_trim_index = sf_start_trim_index(1);
% % % sf_end_trim_index = sf_end_trim_index(1);
% % % 
vf_seq = virtual_force_value.data(vf_start_trim_index+1:vf_end_trim_index,1);
vf_time = virtual_force_value.data(vf_start_trim_index+1:vf_end_trim_index,2);
fx = roundn(virtual_force_value.data(vf_start_trim_index+1:vf_end_trim_index,3),3);
fy = roundn(virtual_force_value.data(vf_start_trim_index+1:vf_end_trim_index,4),3);
fz = roundn(virtual_force_value.data(vf_start_trim_index+1:vf_end_trim_index,5),3);

% % % sf_seq = sensor_force_value.data(sf_start_trim_index+1:sf_end_trim_index,1);
% % % sf_time = sensor_force_value.data(sf_start_trim_index+1:sf_end_trim_index,2);
% % % sfx = roundn(sensor_force_value.data(sf_start_trim_index+1:sf_end_trim_index,3),3);
% % % sfy = roundn(sensor_force_value.data(sf_start_trim_index+1:sf_end_trim_index,4),3);
% % % sfz = roundn(sensor_force_value.data(sf_start_trim_index+1:sf_end_trim_index,5),3);
% % % 

eval_vf_time = vf_time-vf_time(1);

figure ('Name','Virtual Force - Subset')

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
    %ylim([vf_ylim_l,vf_ylim_u])
end

suptitle('Virtual Force  - Subset')


% % % figure ('Name','Sensor Vs Virtual Force - Velocity override 100%')
% % % 
% % % subplot(3,1,1)
% % % plot(vf_time-vf_time(1),fz*-1,'b')
% % % hold on
% % % plot(sf_time-sf_time(1),sfx,'g')
% % % ylabel('sfx/vfz')
% % % legend('-vfz','sfx')
% % % 
% % % subplot(3,1,2)
% % % plot(vf_time-vf_time(1),fy,'b')
% % % hold on
% % % plot(sf_time-sf_time(1),sfy,'g')
% % % ylabel('sfy/vfy')
% % % legend('vfy','sfy')
% % % 
% % % subplot(3,1,3)
% % % plot(vf_time-vf_time(1),fx,'b')
% % % hold on
% % % plot(sf_time-sf_time(1),sfz,'g')
% % % ylabel('sfz/vfx')
% % % legend('vfx','sfz')
% % % 
% % % for i = 1:3
% % %     subplot(3,1,i)
% % %     xlabel('time')
% % % %    ylim([vf_ylim_l,vf_ylim_u])
% % % end
% % % 
% % % suptitle('Sensor Vs Virtual Force - Velocity override 100%')
% % % 
% % % % 20% override 
% % % 
% % % subset_start_time = 100;
% % % subset_end_time = 140;
% % % 
% % % vf_time_trim_1 = eval_vf_time(eval_vf_time>subset_start_time) ;
% % % vf_start_trim_value = vf_time_trim_1(1);
% % % vf_start_trim_index = find( eval_vf_time== vf_start_trim_value);
% % % 
% % % vf_time_trim_2 = eval_vf_time(eval_vf_time<subset_end_time);
% % % vf_end_trim_value = vf_time_trim_2(end);
% % % vf_end_trim_index = find( eval_vf_time== vf_end_trim_value);
% % % 
% % % sf_time_trim_1 = eval_sf_time(eval_sf_time>subset_start_time) ;
% % % sf_start_trim_value = sf_time_trim_1(1);
% % % sf_start_trim_index = find( eval_sf_time == sf_start_trim_value);
% % % 
% % % sf_time_trim_2 = eval_sf_time(eval_sf_time<subset_end_time);
% % % sf_end_trim_value = sf_time_trim_2(end);
% % % sf_end_trim_index = find( eval_sf_time == sf_end_trim_value);
% % % 
% % % vf_start_trim_index = vf_start_trim_index(1);
% % % vf_end_trim_index = vf_end_trim_index(1);
% % % sf_start_trim_index = sf_start_trim_index(1);
% % % sf_end_trim_index = sf_end_trim_index(1);
% % % 
% % % vf_seq = virtual_force_value.data(vf_start_trim_index+1:vf_end_trim_index,1);
% % % vf_time = virtual_force_value.data(vf_start_trim_index+1:vf_end_trim_index,2);
% % % fx = roundn(virtual_force_value.data(vf_start_trim_index+1:vf_end_trim_index,3),3);
% % % fy = roundn(virtual_force_value.data(vf_start_trim_index+1:vf_end_trim_index,4),3);
% % % fz = roundn(virtual_force_value.data(vf_start_trim_index+1:vf_end_trim_index,5),3);
% % % 
% % % sf_seq = sensor_force_value.data(sf_start_trim_index+1:sf_end_trim_index,1);
% % % sf_time = sensor_force_value.data(sf_start_trim_index+1:sf_end_trim_index,2);
% % % sfx = roundn(sensor_force_value.data(sf_start_trim_index+1:sf_end_trim_index,3),3);
% % % sfy = roundn(sensor_force_value.data(sf_start_trim_index+1:sf_end_trim_index,4),3);
% % % sfz = roundn(sensor_force_value.data(sf_start_trim_index+1:sf_end_trim_index,5),3);
% % % 
% % % figure ('Name','Sensor Vs Virtual Force - Velocity override 20%')
% % % 
% % % subplot(3,1,1)
% % % plot(vf_time-vf_time(1),fz*-1,'b')
% % % hold on
% % % plot(sf_time-sf_time(1),sfx,'g')
% % % ylabel('sfx/vfz')
% % % legend('-vfz','sfx')
% % % 
% % % subplot(3,1,2)
% % % plot(vf_time-vf_time(1),fy,'b')
% % % hold on
% % % plot(sf_time-sf_time(1),sfy,'g')
% % % ylabel('sfy/vfy')
% % % legend('vfy','sfy')
% % % 
% % % subplot(3,1,3)
% % % plot(vf_time-vf_time(1),fx,'b')
% % % hold on
% % % plot(sf_time-sf_time(1),sfz,'g')
% % % ylabel('sfz/vfx')
% % % legend('vfx','sfz')
% % % 
% % % for i = 1:3
% % %     subplot(3,1,i)
% % %     xlabel('time')
% % % %    ylim([vf_ylim_l,vf_ylim_u])
% % % end
% % % 
% % % suptitle('Sensor Vs Virtual Force - Velocity override 20%')
% % % 

% 5% override 

% % % subset_start_time = 12;
% % % subset_end_time = 28;
% % % 
% % % vf_time_trim_1 = eval_vf_time(eval_vf_time>subset_start_time) ;
% % % vf_start_trim_value = vf_time_trim_1(1);
% % % vf_start_trim_index = find( eval_vf_time== vf_start_trim_value);
% % % 
% % % vf_time_trim_2 = eval_vf_time(eval_vf_time<subset_end_time);
% % % vf_end_trim_value = vf_time_trim_2(end);
% % % vf_end_trim_index = find( eval_vf_time== vf_end_trim_value);
% % % 
% % % sf_time_trim_1 = eval_sf_time(eval_sf_time>subset_start_time) ;
% % % sf_start_trim_value = sf_time_trim_1(1);
% % % sf_start_trim_index = find( eval_sf_time == sf_start_trim_value);
% % % 
% % % sf_time_trim_2 = eval_sf_time(eval_sf_time<subset_end_time);
% % % sf_end_trim_value = sf_time_trim_2(end);
% % % sf_end_trim_index = find( eval_sf_time == sf_end_trim_value);
% % % 
% % % vf_start_trim_index = vf_start_trim_index(1);
% % % vf_end_trim_index = vf_end_trim_index(1);
% % % sf_start_trim_index = sf_start_trim_index(1);
% % % sf_end_trim_index = sf_end_trim_index(1);
% % % 
% % % vf_seq = virtual_force_value.data(vf_start_trim_index+1:vf_end_trim_index,1);
% % % vf_time = virtual_force_value.data(vf_start_trim_index+1:vf_end_trim_index,2);
% % % fx = roundn(virtual_force_value.data(vf_start_trim_index+1:vf_end_trim_index,3),3);
% % % fy = roundn(virtual_force_value.data(vf_start_trim_index+1:vf_end_trim_index,4),3);
% % % fz = roundn(virtual_force_value.data(vf_start_trim_index+1:vf_end_trim_index,5),3);
% % % 
% % % sf_seq = sensor_force_value.data(sf_start_trim_index+1:sf_end_trim_index,1);
% % % sf_time = sensor_force_value.data(sf_start_trim_index+1:sf_end_trim_index,2);
% % % sfx = roundn(sensor_force_value.data(sf_start_trim_index+1:sf_end_trim_index,3),3);
% % % sfy = roundn(sensor_force_value.data(sf_start_trim_index+1:sf_end_trim_index,4),3);
% % % sfz = roundn(sensor_force_value.data(sf_start_trim_index+1:sf_end_trim_index,5),3);
% % % 
% % % figure ('Name','Sensor Vs Virtual Force - Velocity override 50%')
% % % 
% % % vf_ylim_l = -100; 
% % % vf_ylim_u = 100;
% % % 
% % % subplot(3,1,1)
% % % plot(vf_time-vf_time(1),fz*-1,'b')
% % % hold on
% % % plot(sf_time-sf_time(1),sfx,'g')
% % % ylabel('sfx/vfz')
% % % legend('-vfz','sfx')
% % % 
% % % subplot(3,1,2)
% % % plot(vf_time-vf_time(1),fy,'b')
% % % hold on
% % % plot(sf_time-sf_time(1),sfy,'g')
% % % ylabel('sfy/vfy')
% % % legend('vfy','sfy')
% % % 
% % % subplot(3,1,3)
% % % plot(vf_time-vf_time(1),fx,'b')
% % % hold on
% % % plot(sf_time-sf_time(1),sfz,'g')
% % % ylabel('sfz/vfx')
% % % legend('vfx','sfz')
% % % 
% % % for i = 1:3
% % %     subplot(3,1,i)
% % %     xlabel('time')
% % %     %ylim([vf_ylim_l,vf_ylim_u])
% % % end
% % % 
% % % suptitle('Sensor Vs Virtual Force - Velocity override 50%')
