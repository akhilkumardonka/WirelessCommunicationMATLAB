clc;
clear all;
close all;

% SNR range generation
snr_db=-20:20;
snr_lin=10.^(snr_db./10);

% bit generation
M=10^5;
b_mod=randi([0 1],2,M);

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
    
    H = zeros(M,2);
    H = sqrt(0.5)*[randn(M,2) + 1i*randn(M,2)];
    
    y = zeros(1,M);
    y_tilde = zeros(1,M);
    x_tilde = zeros(2,M);
    
    for m = 1:M
        y(m) = H(m,:)*x(:,m) + sqrt(0.5)*randn(1,1);
        x_tilde(1,m) = -conj(x(2,m));
        x_tilde(2,m) = conj(x(1,m));
        y_tilde(m) = H(m,:)*x_tilde(:,m) + sqrt(0.5)*randn(1,1);
    end
    
    y_p = zeros(2,M);
    for m = 1:M
        y_fin = [y(m);conj(y_tilde(m))];
        big_H = [H(m,1) H(m,2) ; conj(H(m,2)) -conj(H(m,1))];
        w1 = big_H(:,1)/(norm(big_H(:,1)));
        w2 = big_H(:,2)/(norm(big_H(:,2)));
        y_tilde1 = w1'*y_fin;
        y_tilde2 = w2'*y_fin;
        y_p(1,m) = y_tilde1;
        y_p(2,m) = y_tilde2;
    end
    
    x_hat = zeros(2,M);
    for jj=1:M
        if y_p(1,jj)>0
            x_hat(1,jj)=1;
        else
            x_hat(1,jj)=0;
        end
        if y_p(2,jj)>0
            x_hat(2,jj)=1;
        else
            x_hat(2,jj)=0;
        end
        
    end
    error=xor(b_mod,x_hat);
    ber(j)=sum(sum(error))/(2*M);
    
end
semilogy(snr_db,ber);
    
    
    