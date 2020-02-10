clc;
close all;
clear all;

range=-100:0.01:100;
mu=5;
sigma2=300;
r_len=length(range);
% Guassian Theoretical Function
pdf_theo=(1/sqrt(2*pi*sigma2))*exp(-((range-mu*ones(1,r_len)).^2)./(2*sigma2));
plot(range,pdf_theo);
hold on;

samp_size=20000;
data=zeros(1,samp_size);
data=randn(1,samp_size);
data2=sqrt(sigma2)*data+mu.*ones(1,samp_size);

histogram(data2,'Normalization','pdf');