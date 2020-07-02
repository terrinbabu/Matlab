% Kalman filter equation example
clc;
clear all;
close all;

% input

EST = 68; % intial estimate
Eest = 2; % intial error in estimate

Emea = 10; % error in measurement

% measured value (random numbers between 70 and 80)
% generate N random numbers in the interval(a,b) : r = a + (b-a).*rand(N,1)

MEA = 71 + (71-79).*rand(100,1);

for i = 1:100
    KG = Eest / (Eest + Emea);
    EST = EST + KG*(MEA(i) - EST);
    Eest = (1-KG)*Eest;
    kalman_gain(i) = KG;
    estimate(i) = EST;
    error_in_estimate(i) = Eest;
end


figure ('Name','Kalman Gain , Estimate and Error')

subplot(3,1,1)
plot(kalman_gain)
xlabel('time')
ylabel('kalman gain')

subplot(3,1,2)
plot(estimate)
% hold on
% plot(MEA)
xlabel('time')
ylabel('estimate')
% legend('estimate','measure')

subplot(3,1,3)
plot(error_in_estimate)
xlabel('time')
ylabel('error in estimate')

suptitle('Kalman Gain , Estimate and Error')
