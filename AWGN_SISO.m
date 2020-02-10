clc;
clear all;
M=10.^7;
b=rand(1,M);

snr_db=0:15;
snr_lin=10.^(snr_db./10);
ber_theo=zeros(1,length(snr_db));
ber_prac=zeros(1,length(snr_db));

for i=1:length(snr_db)
    ber_theo(i)=qfunc(sqrt(snr_lin(i)));
end

x=zeros(1,M);

for i=1:M
    if b(i)>0.5
        b_mod(i)=1;
    else
        b_mod(i)=0;
    end
end

noise = randn(1,M);

for j=1:length(snr_db)
   p=snr_lin(j);
   for jj=1:M
       if b_mod(jj)==0
           x(jj)=-sqrt(p);
       else
           x(jj)=sqrt(p);
       end
   end
   y=x+noise;
   for jj=1:M
       if y(jj)>0
           x_hat(jj)=1;
       else
           x_hat(jj)=0;
       end
   end
   error=xor(x_hat,b_mod);
   ber_prac(j)=sum(error)/M;
end

semilogy(snr_db,ber_theo);
hold on;
semilogy(snr_db,ber_prac);
title('BER vs log(SNR) - Theoretical and Practical');
ylabel('log(BER)')
xlabel('SNR')
