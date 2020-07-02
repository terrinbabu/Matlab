% FIR filter

function [filtered_signal] = FIR_filter(signal,time,cut_freq,FIR_order)

Fs = 1/mean(diff(time));

Nq = Fs/2;
norm_cut_freq = cut_freq/Nq;
%norm_cut_freq = 2*pi*cut_freq/Fs

[b,a] = fir1(FIR_order,norm_cut_freq);
filtered_signal = filter(b, a, signal);

% % Impulse_respose = fir1(FIR_order,norm_cut_freq);
% % % convolution
% % filtered_signal = conv(signal,Impulse_respose);

end