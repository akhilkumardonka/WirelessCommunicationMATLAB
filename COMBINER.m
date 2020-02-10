clc;
close all;
clear all;

L = [2];

% SNR range generation
snr_db=0:2:20;
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
    
    h_max_abs=max(abs(h).^2);
    %selective combining
    for jj=1:M
        for jjj=1:L(z)
        if(h_max_abs(jj)== abs(h(jjj,jj))^2)
            h_max(jj)=h(jjj,jj);
            n_max(jj)=noise(jjj,jj);
        end
        end
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
       
       for jj=1:M
        y_sec(jj)=h_max(jj)*x(jj)+n_max(jj);
       end
        h_conj=conj(h_max);
   
%   y_hat generation 

   for jj=1:M
       y_hat(jj)=real(h_conj(jj).*y_sec(jj));
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
   ber_prac_sec(j)=sum(error)/M;
   
%    ========================================

       y_tilde_mrc= zeros(1,M);
       w=h./abs(h);
        for jj=1:M
            y(:,jj)=h(:,jj)*x(jj)+noise(:,jj);
            Wopt_mrc=h(:,jj)/norm(h(:,jj));
            Wopt_egc=w(:,jj)/norm(w(:,jj));
            y_tilde_mrc(jj)=Wopt_mrc'*y(:,jj);
            y_tilde_egc(jj)=Wopt_egc'*y(:,jj);
        end
        %Wopt_egc=(h./abs(h))./sqrt(L(z));    
        for jj=1:M
            if real(y_tilde_mrc(jj))>0
                x_hat_mrc(jj)=1;
            else
                x_hat_mrc(jj)=0;
            end
        end
        for jj=1:M
            if real(y_tilde_egc(jj))>0
                x_hat_egc(jj)=1;
            else
                x_hat_egc(jj)=0;
            end
        end
        error_mrc=xor(x_hat_mrc,b_mod);
        ber_prac_mrc(j)=sum(error_mrc)/M;
        error_egc=xor(x_hat_egc,b_mod);
        ber_prac_egc(j)=sum(error_egc)/M;
    end
    semilogy(snr_db,ber_prac_mrc,'LineWidth',2.0);
    hold on;
    semilogy(snr_db,ber_prac_egc,'LineWidth',2.0);
    hold on;
    semilogy(snr_db,ber_prac_sec,'LineWidth',2.0);
    hold on;
end


title('BER vs SNR - Practical (SIMO)');
ylabel('BER');
xlabel('SNR');
legend('Maximum Ratio','EGC','Selective');