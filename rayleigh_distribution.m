clc;
close all;
clear all;

samp_size=1000;
sigma2=1;
x=sqrt(sigma2/2)*randn(1,samp_size);
y=sqrt(sigma2/2)*randn(1,samp_size);
h=x+1j.*y;
a=abs(h);
% a=sqrt(x.^2+y.^2);

histogram(a,'Normalization','pdf');
hold on;

range=0:0.01:100;
r_len=length(range);
pdf_ray_theo=(2.*range).*exp(-range.^2);
plot(range,pdf_ray_theo);