clc;
%close all;
clear all;
N=128;
k=4;
M=10^5;

A = randi([0,1],k,M);
A_mod=A;
for r=1:k
    for c=1:M
        if(A_mod(r,c)==0)
            A_mod(r,c)=-1;
        end
    end
end
a0=A_mod(1,:);
a1=A_mod(2,:);
a2=A_mod(3,:);
a3=A_mod(4,:);

snr_db =-50:5:0;
snr_lin=10.^(snr_db/10);
% H=hadamard(128);
ini1 = [0 0 1 0 1 0 0 1];
ini2 = [0 0 1 0 1 0 0 1];
ini3 = [0 1 1 0 1 0 0 1];
ini4 = [1 0 1 0 1 0 0 1];
pnseq1 = comm.PNSequence('Polynomial','x^8+x^4+x^2+x+1','SamplesPerFrame',128,'InitialConditions',ini1);
pnseq2 = comm.PNSequence('Polynomial','x^8+x^3+x^2+x+1','SamplesPerFrame',128,'InitialConditions',ini2);
pnseq3 = comm.PNSequence('Polynomial','x^8+x^5+x^2+1','SamplesPerFrame',128,'InitialConditions',ini3);
pnseq4 = comm.PNSequence('Polynomial','x^8+x^6+x^3+x+1','SamplesPerFrame',128,'InitialConditions',ini4);
pn1 = pnseq1();
pn2 = pnseq2();
pn3 = pnseq3();
pn4 = pnseq4();
H(1,:) = pn1;
H(2,:) = pn2;
H(3,:) = pn3;
H(4,:) = pn4;


c0=H(1,:);
c1=H(2,:);
c2=H(3,:);
c3=H(4,:);
for i=1:length(snr_db)
    P=snr_lin(i);
    for m=1:M
        x=sqrt(P)*(A_mod(1,m)*c0+A_mod(2,m)*c1+A_mod(3,m)*c2+A_mod(4,m)*c3);
        y=x+randn(1,N);
        d_0=(1/N).*sum(y.*c0);
        if (d_0>0)
            a0_dec(m)=1;
        else
            a0_dec(m)=0;
        end
        
        d_1=(1/N).*sum(y.*c1);
        if (d_1>0)
            a1_dec(m)=1;
        else
            a1_dec(m)=0;
        end
        
        d_2=(1/N).*sum(y.*c2);
        if (d_2>0)
            a2_dec(m)=1;
        else
            a2_dec(m)=0;
        end
        
        d_3=(1/N).*sum(y.*c3);
        if (d_3>0)
            a3_dec(m)=1;
        else
            a3_dec(m)=0;
        end
    end
    error1=xor(A(1,:),a0_dec);
    ber(1,i)=(sum(error1))/M;
    error2=xor(A(2,:),a1_dec);
    ber(2,i)=(sum(error2))/M;
    error3=xor(A(3,:),a2_dec);
    ber(3,i)=(sum(error3))/M;
    error4=xor(A(4,:),a3_dec);
    ber(4,i)=(sum(error4))/M;
end
semilogy(snr_db,ber(1,:))
hold on
semilogy(snr_db,ber(2,:))
hold on
semilogy(snr_db,ber(3,:))
hold on
semilogy(snr_db,ber(4,:))
hold on