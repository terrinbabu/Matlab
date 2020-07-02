close all; clear all; clc;

sine_value = importdata('sine_wave.txt');
fsine_value = importdata('fsine_wave.txt');

sine_time_stamp = sine_value.data(:,1);
sine_value_data = sine_value.data(:,2);

fsine_time_stamp = fsine_value.data(:,1);
fsine_value_data = fsine_value.data(:,2);

% MA Filter
nSamp = 5;
MA_filtered_signal = MA_filter(sine_value_data,nSamp);

% FIR Filter
FIR_order = 32;
cut_freq = 50;
FIR_filtered_signal = LP_FIR_filter(sine_value_data,sine_time_stamp,cut_freq,FIR_order);


% IIR Filter
IIR_order = 10;
cut_freq = 50;
IIR_filtered_signal = LP_butter_filter(sine_value_data,sine_time_stamp,cut_freq,IIR_order);


figure ('Name','Sine waves')
plot(sine_time_stamp,[sine_value_data,MA_filtered_signal,FIR_filtered_signal,IIR_filtered_signal],fsine_time_stamp,fsine_value_data)
xlabel('time') 
ylabel('signal/filtered signal')
legend('req signal','pre filtered signal')

% PSD(sine_value_data,sine_time_stamp);
% PSD(fsine_value_data,fsine_time_stamp);
% PSD(MA_filtered_signal,sine_time_stamp);
% PSD(FIR_filtered_signal,sine_time_stamp);
% PSD(IIR_filtered_signal,sine_time_stamp);
