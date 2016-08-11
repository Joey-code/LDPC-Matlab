% 16QAM 在高斯信道的误比特率计算
function [smld_err_prb] = QAM16_AWGN(SNRindB1)
% close all;
% clear all;

echo on
% SNRindB1=0:2:15;
% SNRindB2=0:0.1:15;
M=16;
k=log2(M);
for i=1:length(SNRindB1)
    %函数开始 
    %smld_err_prb(i)=QAM_awgn_prb(SNRindB1(i));	% simulated error rate
    N=10000;
d=1;				  	% min. distance between symbols
Eav=10*d^2;		 	  	% energy per symbol
snr=10^(SNRindB1(i)/10);	 	  	% SNR per bit (given)
sgma=sqrt(Eav/(8*snr));	  	  	% noise variance
M=16;
% Generation of the data source follows.
for ii=1:N	
  temp=rand;			  	% a uniform R.V. between 0 and 1
  dsource(ii)=1+floor(M*temp);	% a number between 1 and 16, uniform 
end
% Mapping to the signal constellation follows.
mapping=[-3*d 3*d;
	   -d  3*d;
            d  3*d;
	  3*d  3*d;
	 -3*d  d;
	   -d  d;
	    d  d;
	  3*d  d;
 	 -3*d  -d; 
	   -d  -d; 
	    d  -d;
          3*d  -d;
	 -3*d  -3*d;
	   -d  -3*d;
	    d  -3*d;
	  3*d  -3*d];
for ii=1:N
  qam_sig(ii,:)=mapping(dsource(ii),:);
end
% received signal
for ii=1:N
  
  %函数开始 [n(1) n(2)]=gngauss(sgma);
  u=rand;                         	% a uniform random variable in (0,1)       
z=sgma*(sqrt(2*log(1/(1-u))));  	% a Rayleigh distributed random variable
u=rand;                         	% another uniform random variable in (0,1)

n(1) = z*cos(2*pi*u);
n(2) = z*sin(2*pi*u);
%[n(1) n(2)]=[0+z*cos(2*pi*u) 0+z*sin(2*pi*u)];
  %函数结束
  r(ii,:)=qam_sig(ii,:)+n;
end
% detection and error probability calculation
numoferr=0;
for ii=1:N
  % Metric computation follows.
  for jj=1:M
    metrics(jj)=(r(ii,1)-mapping(jj,1))^2+(r(ii,2)-mapping(jj,2))^2;
  end
  [min_metric decis] = min(metrics);
  if (decis~=dsource(ii)),
    numoferr=numoferr+1;
  end
end
smld_err_prb(i)=numoferr/(N);	

    %函数结束
 
  echo off;
end
% semilogy(SNRindB1,smld_err_prb,'*');
% hold on
% echo on ;
%for i=1:length(SNRindB2)
  %SNR=exp(SNRindB2(i)*log(10)/10);    	% signal-to-noise ratio
  % theoretical symbol error rate
  %函数开始 theo_err_prb(i)=4*Qfunct(sqrt(3*k*SNR/(M-1)));  
  %theo_err_prb(i)=4*(1/2)*erfc(sqrt(3*k*SNR/(M-1))/sqrt(2));
  %函数结束
  
  %echo off ;
%end
%echo on ;
% Plotting commands follow.

% hold on
% semilogy(SNRindB2,theo_err_prb);