clc;
clear all;
tic;
m=128;n=256;
R=(n-m)/n
frame=500;
N=frame*n;
Eb_N0=[0:0.5:10];

H=getH(m,n);
[G,valid]=H2G(H); 
while valid==0                  
H = getH(m,n);               
[G,valid]=H2G(H); 
end


  for i=1:length(Eb_N0)
    %sigma_2=1/(2*(10^(Eb_N0(i)/10))*R); 
    sigma=sqrt(1./(2*10^(Eb_N0(i)/10)*R));
    %N0 = 1/(exp(Eb_N0(i)*log(10)/10));
    %sigma=sqrt(N0/2)
    ber0(i)=0;
    %sigma_2
    %i
    for num=1:frame  
        num
        %x=round(rand(1,n-m));
        x = (sign(randn(1,size(G,1)))+1)/2; % random bits
        y = mod(x*G,2);                     % coding 
        bpskmod = 2*y-1;                          %BPSK modulation
        h = 1/sqrt(2)*[randn(1,n) + j*randn(1,n)]; % Rayleigh channel
        z =  h.*y + sigma*randn(size(bpskmod));%����˹���������Ϣ���У�Ҳ��׼�����������
        yHat = z./h;   % equalization
        ipHat = real(yHat)>0;
        %z=bpskmod + sigma*randn(1,size(G,2));   % AWGN transmission
        f1=1./(1+exp(-2* ipHat/(sigma^2)));         % likelihoods
        f0=1-f1;
        [z_hat, success, k] = ldpc_decode(ipHat,f0,f1,H);
        x_hat = z_hat(size(G,2)+1-size(G,1):size(G,2));
        x_hat = x_hat';
        err_max=find(x~=x_hat);             %Ѱ�Ҵ�����Ϣλ
        num_eer=length(err_max)           %���������Ϣλλ��
        ber(i)=num_eer/n;               %�������������BER
        ber0(i)=ber(i)+ber0(i);
        %ber(i);
        %ber0(i);
        %k;
        
    end %for num
    ber0(i)=ber0(i)/frame
    
end %for i
semilogy(Eb_N0,ber0,'b-o');
