clc;
clear all;
m=256;n=512;
R=(n-m)/n
frame=1000;
Eb_N0=[0:1:10];

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
        %x= rand(1,n)>0.5;
        %x=round(rand(1,m));
        x = (sign(randn(1,size(G,1)))+1)/2; % random bits
        y = mod(x*G,2);                     % coding 
        bpskmod = 2*y-1;                          %BPSK modulation
        h = 1/sqrt(2)*abs([randn(1,n) + j*randn(1,n)]); % Rayleigh channel
        z =h.*bpskmod + sigma*randn(size(bpskmod));%��Rayleigh�ŵ�����Ϣ���У�Ҳ��׼�����������
        %z=bpskmod + sigma*randn(1,size(G,2));   % AWGN transmission
        f1=1./(1+exp(-2*z*0.8862/(sigma^2)));         % likelihoods
        f0=1-f1;
        [z_hat, success, k] = ldpc_decode(z,f0,f1,H);
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
hold on
%%%%%%%%%%%%%%%%%%%% QPSK %%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
N = 1000;
iter = 100;
expSymErrs = 10;
EbNo = 0:1:10;


for t = 1:length(EbNo)
    idx = 1;
    errSym = 0;
    while ((idx <= iter) || (sum(errSym) <= expSymErrs))
        s = 1/sqrt(2)*(randsrc(N,1) + j*randsrc(N,1));
        % Compute EsNo
        EsNo = EbNo(t) + 10*log10(2);
        EsNolin = 10.^(EsNo/10);
        No = 1/EsNolin; % Assuming Es = 1
        var = No;     %
        n = sqrt(var/2)*(randn(N,1) + j*randn(N,1));
        x = s + n;
        % Detection follows
        alph = [1+j, 1-j, -1+j, -1-j]/sqrt(2);
        detx = zeros(N,1);
        for tt = 1:N
            [val,id] = min(abs(x(tt) - alph));
            detx(tt) = alph(id);
        end
        errSym(idx) = length(find(s ~= detx));
        idx = idx + 1;
    end
    pe(t) = mean(errSym)/N;
end
semilogy(EbNo,pe,'k-d')
hold on;

%%%%%%%%%%%%%%%%%%%16QAM %%%%%%%%%%%%%
 clc;
 SNRindB1=0:1:10;
 
 
 %***��ÿ����������ȵ�ͨ��ϵͳ������BER��ֵ
 [smld_err_prb] = QAM16_AWGN(SNRindB1);
 pb1 = smld_err_prb;

 
semilogy(SNRindB1, pb1,'r-*');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlabel('SNR(dB)');
ylabel('BER');
legend('BPSK','QPSK','16QAM',1);
grid on;
hold off;
        
        
        
        
        
        
        
        