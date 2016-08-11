function modout=modulation(mod_in,mod_mode)
%%*************************************************************************
%%Function information:
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
%% 
%%-------------------------------------------------------------------------

%% Function discription:
%%-------------------------------------------------------------------------
%%��������ĵ��Ʒ�ʽ������������MOD_IN���е��ƣ��ֱ����BPSK, QPSK, !6QAM, 64QAM��
%%��ɶ�����ͼ��ӳ�䣬���ΪY.ת���ķ���Ϊ����д��ʮ��������´�0 ��N-1
%%��NΪ����ͼ�ĵ���������Ӧ���������ꣻ�ٽ�����Ķ���������ת��Ϊ��Ӧ��
%%ʮ���ƣ��Բ��ķ��������Ӧ��ĸ������꣬��Ϊ����ӳ���Ľ����

%%The OFDM subcarriers shall be modulated by using BPSK, QPSK, 16-QAM, or 64-QAM modulation,
%%depending on the RATE requested. The encoded and interleaved binary serial input data shall
%%be divided into groups of N BPSC (1, 2, 4, or 6) bits and converted into complex numbers 
%%representing BPSK, QPSK, 16-QAM, or 64-QAM constellation points. The conversion shall be
%%performed according to Gray-coded constellation mappings, with the input bit MOD_IN. The
%%output values,MOD_OUT are formed by multiplying the resulting (I+jQ) value by a normalization 
%%factor K MOD , as described in equation d = (I + jQ) �� K MOD ) The normalization factor,K MOD,
%%depends on the base modulation mode. Note that the modulation type can be different from the
%%start to the end of the transmission, as the signal changes from SIGNAL to DATA. The purpose
%%of the normalization factor is to achieve the same  average power for all mappings.In practical
%%implementations, an approximate value of the normalization factor can be used, as long as 
%%the device conforms with the modulation accuracy requirements .

%%-------------------------------------------------------------------------

%% Input: 
%%-------------------------------------------------------------------------
%%  mod_in:����Ķ���������(The sequence to be modulated)
%%-------------------------------------------------------------------------
%% Output:
%%-------------------------------------------------------------------------
%% mod_out:����ͼӳ���õ��ĵ��Ƹ������(The output after modulation)
%%-------------------------------------------------------------------------
%% Global Variable:
%%  g_RT (the vector which contains the modulation mode)
%%-------------------------------------------------------------------------
%% Z :ѡ����Ʒ�ʽ�Ĳ��� (the parameter to choose the modulation mode)
%%
%% R :��������������������У���һ��Ҫ�󣩺�Ľ��,���磺��16QAM��Ҫ���������е���Ϊ
%% 4�У�length(g_MOD_IN_16QAM )/4  �еľ���(Reshape the input binary sequence to be
%% matrix of n-row,m-column .For example,16QAM,will be reshaped into 4-row,
%%length(g_MOD_IN_16QAM))/4 column )
%%
%% B2D :��������ʮ����ת����Ľ�� (convert the binary sequence to  dec )
%%
%% Temp:����ͼ���� (the  constellation)
%%-------------------------------------------------------------------------
%%********************************************************
%system_parameters  
switch (mod_mode)
case 2 
    for i=1:length(mod_in)
       if  mod_in(i)==0
           mod_out(i)=-1;
       else mod_out(i)=1;%�������ͼ��ӳ�� (mapping)
       end
    end
%disp(mod_out)  %������

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 4
    mod_out=zeros(1,length(mod_in)/2);
    R=reshape(mod_in,2,length(mod_in)/2); %����������ת��Ϊ(2,length(x)/2)�ľ���
    B2D=bi2de(R','left-msb')+1;  %��������תΪʮ����,ע���1����Ϊmatlabû��a(0)����Ǵ�a(1)��ʼ
    Temp=[-1-j  -1+j  1-j   1+j];
    for i=1:length(mod_in)/2
       mod_out(i)=Temp(B2D(i))/sqrt(2);%��һ��
       modout(i)=real( mod_out(i));
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 16
    mod_out=zeros(1,length(mod_in)/4);
    R=reshape(mod_in,4,length(mod_in)/4);
    B2D=bi2de(R','left-msb')+1;
    Temp=[-3-3*j   -3-j   -3+3*j   -3+j ...   
          -1-3*j   -1-j   -1+3*j   -1+j ...
           3-3*j    3-j    3+3*j    3+j ...
           1-3*j    1-j    1+3*j    1+j];
   for i=1:length(mod_in)/4
       mod_out(i)=Temp(B2D(i))/sqrt(10);
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 64
   mod_out=zeros(1,length(mod_in)/6);
   R=reshape(mod_in,6,length(mod_in)/6);
   B2D=bi2de(R','left-msb')+1;
   Temp=[-7-7*j  -7-5*j  -7-j  -7-3*j  -7+7*j  -7+5*j  -7+j  -7+3*j...  
         -5-7*j  -5-5*j  -5-j  -5-3*j  -5+7*j  -5+5*j  -5+j  -5+3*j... 
         -1-7*j  -1-5*j  -1-j  -1-3*j  -1+7*j  -1+5*j  -1+j  -1+3*j...
         -3-7*j  -3-5*j  -3-j  -3-3*j  -3+7*j  -3+5*j  -3+j  -3+3*j...
          7-7*j   7-5*j   7-j   7-3*j   7+7*j   7+5*j   7+j   7+3*j...
          5-7*j   5-5*j   5-j   5-3*j   5+7*j   5+5*j   5+j   5+3*j...
          1-7*j   1-5*j   1-j   1-3*j   1+7*j   1+5*j   1+j   1+3*j...
          3-7*j   3-5*j   3-j   3-3*j   3+7*j   3+5*j   3+j   3+3*j ];
   for i=1:length(mod_in)/6
       mod_out(i)=Temp(B2D(i))/sqrt(42);
   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
otherwise
    disp('Error! Please input again');
end
   