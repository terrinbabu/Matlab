% MA filter

function [filtered_signal] = MA_filter(signal,n)
coeffMA = ones(1,n)/n;
filtered_signal = filter(coeffMA, 1, signal);
end