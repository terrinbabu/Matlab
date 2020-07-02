% LP butter IIR filter

function [filtered_signal] = LP_butter_filter(signal,time,cut_freq,IIR_order)

Fs = 1/mean(diff(time));

Nq = Fs/2;
norm_cut_freq = cut_freq/Nq;
%norm_cut_freq = 2*pi*cut_freq/Fs;

[b,a] = butter(IIR_order,norm_cut_freq);
filtered_signal = filter(b, a, signal);

% % Impulse_respose = butter(IIR_order,norm_cut_freq);
% % % convolution
% % filtered_signal = conv(signal,Impulse_respose);