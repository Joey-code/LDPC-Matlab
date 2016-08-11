function [n1,n2]=gnray(m,sgma)
% [n1 n2]=gnray(m,sgma)
% [n1 n2]=gnray(sgma)
% [n1 n2]=gnray
% GNRAY generates two independent rayleigh variables with 
% mean m & S.D sgma
% if one of the input argement missing it takes mean as zero
% if neither mean nor variance is given, it generates two 
%standard rayleigh random variables

if nargin==0
    m=0;
    sgma=1;
elseif nargin==1
    sgma=m;  
    m=0;
end
% rayleigh function
u=2*pi*rand(20,1);
n=m+(sqrt(sgma/20))*randn(20,1);
nr1=sum(n.*exp(j*u));
nr2=abs(nr1);
nr3=angle(nr1);
n1=nr2.*cos(nr3);
n2=nr2.*sin(nr3);