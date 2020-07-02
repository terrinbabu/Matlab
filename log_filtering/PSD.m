% Power Spectral Density

function [] = PSD(signal,time)

fsamp = 1/mean(diff(time));
Nfft = 2.^(fix(log2(length(signal))));
[Pxx,f] = pwelch(signal,gausswin(Nfft),Nfft/2,Nfft,fsamp);

figure ('Name', strcat('PSD_',inputname(1)))
plot(f,Pxx);
ylabel('PSD'); xlabel('Frequency (Hz)');
grid on;

[~,loc] = max(Pxx);
FREQ_ESTIMATE = f(loc);
title(['Frequency estimate = ',num2str(FREQ_ESTIMATE),' Hz']);
end