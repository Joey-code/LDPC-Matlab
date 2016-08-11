function  demod_out = demodulation(demod_in,mod_mode)
%%*************************************************************************
%%Function information:
%%���ŵ��źŵĸ���ֵ�����Ӧ�Ķ���������(In this program,the channel value are demodulated.)
%%-------------------------------------------------------------------------
%%First time   : 3/25/2002                                           
%%Newest modified time:6/20/2002                                         
%%Programmer:Xuewei Mao
%%Version: 0.2
%%-------------------------------------------------------------------------
%%*************************************************************************

%%*************************************************************************
%% Reference:
%%-------------------------------------------------------------------------
%% 
%%-------------------------------------------------------------------------

%% Note:
%%-------------------------------------------------------------------------
%% Here,the method for the four modulation modes is the same ,so see the note
%%for QPSK, and the note for the other three is omitted. 
%%-------------------------------------------------------------------------

%% Function discription:
%%-------------------------------------------------------------------------
%% ���ŵ��źŵĸ���ֵ�ҳ���Ӧ�Ķ��������С��������£����ŵ�ֵ��
%%�����ֵ������ͼ�����е�ľ��룬�ҳ�������С�ĵ㣬�õ��Ӧ��
%%���������м�Ϊ���ŵ���Ӧ�Ľ�������
%% In this program,the channel value are demodulated.
%% The process: firstly,calculate  the distance between the channel value
%% and all the constellation points.then,find the shortest distance,thus ,the point
%% is the one.At last,put out the binary sequence symbolizing the point.   
%%-------------------------------------------------------------------------

%% Input: 
%%-------------------------------------------------------------------------
%% demod_in :�����ŵ�ֵ  (the input channel value)
%%-------------------------------------------------------------------------
%% Output:
%%-------------------------------------------------------------------------
%% demod_out:������ (demodulation output(binary) )
%%-------------------------------------------------------------------------
%% Global Variable:
%%-------------------------------------------------------------------------
%%  g_RT
%%-------------------------------------------------------------------------
%%*************************************************************************
%system_parameters  
switch (mod_mode)
case 2
    %BPSK modulation
    demod_out=zeros(1,length(demod_in));
     for i=1:length(demod_in)
       if abs( demod_in(i)-1)-abs( demod_in(i)+1)<=0
          demod_out(i)=1;
       else 
          demod_out(i)=0;
       end
     end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 4
    %QPSK modulation
d=zeros(4,length(demod_in));    % 'd ' is the distance between the channel value and the constellation point  
m=zeros(1,length(demod_in));
temp=[-1-j  -1+j  1-j   1+j]/sqrt(2);
for i=1:length(demod_in)
    for n=1:4
        d(n,i)=(abs(demod_in(i)*sqrt(2)-temp(n))).^2;%���ŵ�ֵ�������ֵ������ͼ�����е�ľ���
    end                %calculate the distance between the channel and the constellation points
    [min_distance,constellation_point] = min(d(:,i)) ; %D(:,i)=sort(d(:,i));%���� (sort the distance in  ascending order.)
    m(i) = constellation_point;
end 
A=de2bi([0:3],'left-msb');%д��0��N-1(NΪ����ͼ����) 
%(converts a nonnegative  decimal vector [0,3] to a binary  matrix A.)

for i=1:length(demod_in)
    DEMOD_OUT(i,:)=A(m(i),:);
   % ��Сֵ��Ӧ����������ŵĶ����Ʊ�ʾ��Ϊ������
end
demod_out=reshape(DEMOD_OUT',1,length(demod_in)*2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 16
d=zeros(16,length(demod_in));
b=ones(1,length(demod_in));
m=zeros(1,length(demod_in));
temp=[-3-3*j   -3-j   -3+3*j   -3+j ...   
      -1-3*j   -1-j   -1+3*j   -1+j ...
       3-3*j    3-j    3+3*j    3+j ...
       1-3*j    1-j    1+3*j    1+j]/sqrt(10);
for i=1:length(demod_in)
    for n=1:16
        d(n,i)=(abs(demod_in(i)-temp(n))).^2;
    end 
    [min_distance,constellation_point] = min(d(:,i)) ; %D(:,i)=sort(d(:,i));%���� (sort the distance in  ascending order.)
end 

A=de2bi([0:15],'left-msb');
for i=1:length((demod_in))
    DEMOD_OUT(i,:)=A(m(i),:);
end
demod_out=reshape(DEMOD_OUT',1,length(demod_in)*4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 64
    %64_QAM
d=zeros(64,length(demod_in));
D=zeros(64,length(demod_in));
m=zeros(1,length(demod_in));
temp=[-7-7*j  -7-5*j  -7-j  -7-3*j  -7+7*j  -7+5*j  -7+j  -7+3*j...  
      -5-7*j  -5-5*j  -5-j  -5-3*j  -5+7*j  -5+5*j  -5+j  -5+3*j... 
      -1-7*j  -1-5*j  -1-j  -1-3*j  -1+7*j  -1+5*j  -1+j  -1+3*j...
      -3-7*j  -3-5*j  -3-j  -3-3*j  -3+7*j  -3+5*j  -3+j  -3+3*j...
       7-7*j   7-5*j   7-j   7-3*j   7+7*j   7+5*j   7+j   7+3*j...
       5-7*j   5-5*j   5-j   5-3*j   5+7*j   5+5*j   5+j   5+3*j...
       1-7*j   1-5*j   1-j   1-3*j   1+7*j   1+5*j   1+j   1+3*j...
       3-7*j   3-5*j   3-j   3-3*j   3+7*j   3+5*j   3+j   3+3*j ]/sqrt(42);
for i=1:length(demod_in)
    for n=1:64
        d(n,i)=(abs(demod_in(i)*sqrt(42)-temp(n))).^2;
    end 
         [min_distance,constellation_point] = min(d(:,i)) ; %D(:,i)=sort(d(:,i));%���� (sort the distance in  ascending order.)
end 
A=de2bi([0:63],'left-msb');
for i=1:length(demod_in)
    DEMOD_OUT(i,:)=A(m(i),:);
end
demod_out=reshape(DEMOD_OUT',1,length(demod_in)*6);
otherwise
    disp('Error! Please input again');        
end
      
