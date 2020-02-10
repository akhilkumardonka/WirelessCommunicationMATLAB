clc;
clear all;
close all;


% SNR range generation
snr_db=0:30;
snr_lin=10.^(snr_db./10);

% bit generation
M=10.^5;
b=randn(1,M);
b_mod=zeros(1,length(b));

for i=1:M
    if b(i)>0.5
        b_mod(i)=1;
    else
        b_mod(i)=0;
    end
end

% noise generation
noise = randn(1,M);

% BER_prac and theo initialization
ber_prac=zeros(1,length(snr_db));
ber_theo=zeros(1,length(snr_db));

% received vector initialization
y=zeros(1,M);

% rayleigh fading channel generation
sigma2=1;
x1=sqrt(sigma2/2)*randn(1,M);
y1=sqrt(sigma2/2)*randn(1,M);
h=x1+1j.*y1;

x=zeros(1,M);

b_hat=zeros(1,M);

for j=1:length(snr_db)
   p=snr_lin(j);
   for jj=1:M
       if b_mod(jj)==0
           x(jj)=-sqrt(p);
       else
           x(jj)=sqrt(p);
       end
   end
   
   for jj=1:M
       y(jj)=h(jj).*x(jj)+noise(jj);
   end
   h_conj=conj(h);
   
%   y_hat generation 

   for jj=1:M
       y_hat(jj)=real(h_conj(jj).*y(jj));
   end
   
   % b_hat generation
   for jj=1:M
       if y_hat(jj)>0
           x_hat(jj)=1;
           b_hat(jj)=1;
       else
           x_hat(jj)=-1;
           b_hat(jj)=0;
       end
   end
   
   % ber_prac generation
   error=xor(b_hat,b_mod);
   ber_prac(j)=sum(error)/M;
end

% ber_theo generation
for i=1:length(snr_db)
    ber_theo(i)=0.5*(1- sqrt(2*snr_lin(i)./(2*snr_lin(i)+2)));
end

% plotting
semilogy(snr_db,ber_theo);
hold on;
semilogy(snr_db,ber_prac);
title('BER vs log(SNR) - Theoretical and Practical(SISO)');
ylabel('log(BER)')
xlabel('SNR')