close all; clc;

velocity = 0.001:0.01:1;
columb = 13.80;
viscous_1o = 21.48;
viscous_2o = -8.504;

for i=1:length(velocity)
    fric_force(i) = columb +  viscous_1o*velocity(i) + viscous_2o*(velocity(i)^2);
end
plot(velocity,fric_force)