clc;
clear all;
close all;

% SNR range generation
snr_db=-20:5:20;
snr_lin=10.^(snr_db./10);

t=2;
r=2;
% bit generation
M=10^5;
b=rand(2,M);

b_mod=zeros(2,M);

for ii=1:2
for i=1:M
    if b(ii,i)>0.5
        b_mod(ii,i)=1;
    else
        b_mod(ii,i)=0;
    end
end
end

for j=1:length(snr_db)
    p=snr_lin(j);
    for k=1:2
        for jj=1:M
            if b_mod(k,jj)==0
                x(k,jj)=-sqrt(p);
            else
                x(k,jj)=sqrt(p);
            end
        end
    end
    H = zeros(r,t,M);
    H = sqrt(0.5)*[randn(r,t,M) + 1i*randn(r,t,M)];
    y = zeros(2,M);
    
    for m = 1:M
        y(:,m) = H(:,:,m)*x(:,m) + sqrt(0.5)*randn(2,1);
        x_zf(:,m)=inv(H(:,:,m)'*H(:,:,m))*H(:,:,m)'*y(:,m);
        for mm=1:2
            if real(x_zf(mm,m))>0
                x_map(mm,m)=1;
            else
                x_map(mm,m)=0;
            end
        end
    end
    
    error=xor(b_mod,x_map);
    ber(j)=sum(sum(error))/(2*M);
end

semilogy(snr_db,ber);