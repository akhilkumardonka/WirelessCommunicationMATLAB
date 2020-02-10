clc;
close all;
clear all;

k=4;
N=128;
M=10^5;

% SNR range generation
snr_db=-20:2:0;
snr_lin=10.^(snr_db./10);
ber=zeros(k,length(snr_db));

% bit generation
b=zeros(k,M);
b_mod=zeros(k,M);

for i= 1:k
        b(i,:)=rand(1,M);
        for ii=1:M
            if b(i,ii)>0.5
                b_mod(i,ii)=1;
            else
                b_mod(i,ii)=0;
            end
        end
end
x=zeros(k,M);

%Walsh code
h=hadamard(N);
c= zeros(k,N);
for i=1:k
    c(i,:)=h(i,:);
end

% noise generation
noise= randn(1,N);

for jj=1:k
       for jjj=1:M
           if b_mod(jj,jjj)==0
               x(jj,jjj)=-1;
           else
               x(jj,jjj)=1;
           end
       end
   end

for j=1:length(snr_db)
   p=snr_lin(j);
   

       for jj=1:M
            a=x(1,jj).*c(1,:)+x(2,jj).*c(2,:)+x(3,jj).*c(3,:)+x(4,jj).*c(4,:);
            y=sqrt(p).*a+noise;
            
            for jjj=1:k
                 d(jjj)=(sum(y.*c(jjj,:)))/N;
                 if(d(jjj)>0)
                     b_hat(jjj,jj)=1;
                 else
                     b_hat(jjj,jj)=0;
                 end
            end
       end
       

       
       for jj=1:k
            error(jj,:)=xor(b_hat(jj,:),b_mod(jj,:));
            ber(jj,j)=sum(error(jj,:))/M;
       end
end

semilogy(snr_db,ber(1,:));

