% BS Cheby IIR filter

function [filtered_signal] = BS_cheby_filter(signal,time,cut_freq1,cut_freq2,filter_order,ripple_db,freq_graph)

Fs = 1/mean(diff(time));

Nq = Fs/2;
norm_cut_freq1 = cut_freq1/Nq;
norm_cut_freq2 = cut_freq2/Nq;

[b,a] = cheby1(filter_order,ripple_db,[ norm_cut_freq1,norm_cut_freq2 ],'stop');
filtered_signal = filter(b, a, signal);

if (freq_graph ~= 0)
    figure ('Name','BS_cheby_filter freq_response')
    freqz(b,a)
end