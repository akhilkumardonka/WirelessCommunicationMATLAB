clc;
close all;
clear all;

L = [2,3,4];

% SNR range generation
snr_db=0:30;
snr_lin=10.^(snr_db./10);

% bit generation
M=10^5;
b=rand(1,M);
b_mod=zeros(1,length(b));

for i=1:M
    if b(i)>0.5
        b_mod(i)=1;
    else
        b_mod(i)=0;
    end
end

for z=1:length(L)
    
    % noise generation
    noise = zeros(L(z),M);
    for i=1:L(z)
        noise(i,:) = sqrt(0.5)*(randn(1,M) + 1i*randn(1,M));
    end

    % BER_prac and theo initialization
    ber_prac=zeros(1,length(snr_db));
    ber_theo=zeros(1,length(snr_db));

    % received vector initialization
    y=zeros(L(z),M);

    % rayleigh fading channel generation
    sigma2=1;
    for i=1:L(z)
        x1(i,:)=sqrt(sigma2/2)*randn(1,M);
        y1(i,:)=sqrt(sigma2/2)*randn(1,M);
        h(i,:)=x1(i,:)+1j.*y1(i,:);
    end

    x=zeros(1,M);

    for j=1:length(snr_db)
       p=snr_lin(j);
       for jj=1:M
           if b_mod(jj)==0
               x(jj)=-sqrt(p);
           else
               x(jj)=sqrt(p);
           end
       end

       y_tilde= zeros(1,M);
        for jj=1:M
            y(:,jj)=h(:,jj)*x(jj)+noise(:,jj);
            Wopt=h(:,jj)/norm(h(:,jj));
            y_tilde(jj)=Wopt'*y(:,jj);
        end
        for jj=1:M
            if real(y_tilde(jj))>0
                x_hat(jj)=1;
            else
                x_hat(jj)=0;
            end
        end
        error=xor(x_hat,b_mod);
        ber_prac(j)=sum(error)/M;
    end
    semilogy(snr_db,ber_prac,'LineWidth',2.0);
    hold on;
end

for j=1:length(L)
    ber = zeros(0,length(snr_db));
    for i=1:length(snr_db)
        lamda= sqrt(2*snr_lin(i)/(2*snr_lin(i)+2));
        sum=0;
            for l=0:L(j)
                sum=sum+nchoosek(L(j)+l-1,l)*((1+lamda)/2)^l;
            end
        ber(i) =(((1-lamda)/2)^L(j))*sum;
    end
semilogy(snr_db,ber,'LineWidth',2.0);
hold on
end

title('BER vs SNR - Theoretical and Practical (SIMO)');
ylabel('BER');
xlabel('SNR');
legend('Analytical L=2','Analytical L=3','Analytical L=4','Theoretical L=2','Theoretical L=3','Theoretical L=4');