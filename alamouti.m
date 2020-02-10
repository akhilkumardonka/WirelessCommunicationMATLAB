clc;
close all;

M=10.^4;
snr_db=-20:5:20;
snr_lin=10.^(snr_db./10);

sigma2=1;

x=randn(2,2,M);
yl=randn(2,2,M);
h=sqrt(sigma2/2).*(x+1j.*yl);

noise=sqrt(sigma2/2).*(randn(2,M)+1j.*randn(2,M));

b_mod=randi([0 1],2,M);
data=zeros(2,M);
x_hat=zeros(2,M);
b_hat=zeros(2,M);
ber_prac=zeros(1,length(snr_db));

for j=1:length(snr_db)
    p=snr_lin(j);
    for jj=1:M
        for jjj=1:2
            if b_mod(jjj,jj)==1
                x_hat(jjj,jj)=sqrt(p);
            else
                x_hat(jjj,jj)=-sqrt(p);
            end
        end
        
        h_tildae=[h(1,:,jj); h(2,:,jj); conj(h(1,2,jj)) -conj(h(1,1,jj));conj(h(2,2,jj)) -conj(h(2,1,jj))];
        noise_tildae=[noise(1,jj);noise(2,jj);conj(noise(1,jj));conj(noise(2,jj))];
        
        y=h_tildae*x_hat(:,jj)+noise_tildae;
        
        h_one=h_tildae(:,1);
        h_two=h_tildae(:,2);
        
        b_hat(1)=(h_one)'*y;
        b_hat(2)=(h_two)'*y;
        
        for jjj=1:2
            if real(b_hat(jjj))>0
                data(jjj,jj)=1;
            else
                data(jjj,jj)=0;
            end
        end
        
        final=xor(data,b_mod);
        error=sum(sum(final))/(2*M);
        
    end
    ber_prac(j)=error;
end

semilogy(snr_db,ber_prac);


