clear all;
clc;

v = 16.67;
fc = 2*10.^9;
c=3*10.^8;

fdmax=(v/c)*fc;
disp('fdmax ');disp(fdmax);
tc = 1/(4*fdmax);
k=-100:1:100;
fs=1/(1*tc);
df=fs/length(k);
f=-fs/2:df:fs/2-df;

arg=(pi/2)*k;

J = besselj(0,arg);
f_j=fft(J);
% plot(k,J,'LineWidth',1.5)
% plot(k,abs(f_j),'LineWidth',1.5)
% plot(k,fftshift(abs(f_j)),'LineWidth',1.5)
plot(f,abs(f_j));
title('Jakes Correlation as function of normalised form')
xlabel('K')
ylabel('J_0(arg)')

