% % Linear phase filters - preserve shape of a filtered signal
% % https://www.youtube.com/watch?v=xPTe7ZWLVhQ
% % https://dadorran.wordpress.com/2014/10/01/linear-phase-filters-why-they-are-used/


close all ; clc;

%% synthesise a signal

% time vector 

fs = 100; % sampling frequency
T = 1/fs; %sampling interval
N = 2000; %length of signal being synthesised
n = 0:N-1; %samples of the signal
t = n*T; % time vector

% signal

A1 = 1;
f1 = 10;
phase1 = 0;

A2 = 0.5;
f2 = 20;
phase2 = 1.4;

x = A1*cos(2*pi*f1*t + phase1) + 0.5*cos(2*pi*f2*t + phase2);



% noise

ns = randn(1,length(x))*3;

%filter the noise to synthesise band limited noise
butter_filter_order = 5;
cut_freq1 = 14;
cut_freq2 = 16;
ns_filtered = BP_butter_filter(ns,t,cut_freq1,cut_freq2,butter_filter_order,0);

x_ns = x + ns_filtered;
  

%% using an IIR filter (non-linear phase)

cheby_filter_order = 10;
ripple_db = 0.5;
cut_freq1 = 12;
cut_freq2 = 18;

y_iir = BS_cheby_filter(x_ns,t,cut_freq1,cut_freq2,cheby_filter_order,ripple_db,1);

%% using an FIR filter (linear phase)

fir_filter_order = 100;
cut_freq1 = 12;
cut_freq2 = 18;

y_fir = BS_fir_filter(x_ns,t,cut_freq1,cut_freq2,fir_filter_order,1);

%% Plots

PSD(x,t)
PSD(x_ns,t)
PSD(y_iir,t)
PSD(y_fir,t)

% % plot_range = N/2-100:N/2+100;
% % 
% % figure ('Name','Synthesised Signals')
% % plot(t(plot_range),x(plot_range),'b');
% % hold on;
% % noisy_x = plot(t(plot_range),x_ns(plot_range),'r');
% % hold on;
% % iir_filtered_x = plot(t(plot_range),y_iir(plot_range),'g');
% % hold on
% % fir_filtered_x = plot(t(plot_range),y_fir(plot_range),'k');
% % xlabel('Time (seconds)')
% % ylabel('Amplitude')
% % title('Synthesised Signals')
% % legend('clean signal', 'noisy signal','iir filtered','fir filtered')
% % axis tight
% % pause 
% % set(noisy_x,'visible', 'off')
% % pause
% % set(iir_filtered_x,'visible', 'off')
