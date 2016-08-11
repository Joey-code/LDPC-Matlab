function [x_hat, success, k] = ldpc_decode(y,f0,f1,H)
% decoding of binary LDPC as in Elec. Letters by MacKay&Neal 13March1997
% For notations see the same reference.
% function [x_hat, success, k] = ldpc_decode(y,f0,f1,H)
% outputs the estimate x_hat of the ENCODED sequence for
% the received vector y with channel likelihoods of '0' and '1's
% in f0 and f1 and parity check matrix H. Success==1 signals
% successful decoding. Maximum number of iterations is set to 100.
% k returns number of iterations until convergence.
%
% Example:
% We assume G is systematic G=[A|I] and, obviously, mod(G*H',2)=0 
%         sigma = 1;                          % AWGN noise deviation
%         x = (sign(randn(1,size(G,1)))+1)/2; % random bits
%         y = mod(x*G,2);                     % coding 
%         z = 2*y-1;                          % BPSK modulation
%         z=z + sigma*randn(1,size(G,2));     % AWGN transmission
%
%         f1=1./(1+exp(-2*z/sigma^2));        % likelihoods
%         f0=1-f1;
%         [z_hat, success, k] = ldpc_decode(z,f0,f1,H);
%         x_hat = z_hat(size(G,2)+1-size(G,1):size(G,2));
%         x_hat = x_hat'; 

%   Copyright (c) 1999 by Igor Kozintsev igor@ifp.uiuc.edu
%   $Revision: 1.1 $  $Date: 1999/07/11 $
%   fixed high-SNR decoding

[m,n] = size(H); 
if m>n, H=H'; 
 [m,n] = size(H); 
end
if ~issparse(H) % make H sparse if it is not sparse yet如果不是稀疏矩阵，将H变为稀疏矩阵
   [ii,jj,sH] = find(H);
   H = sparse(ii,jj,sH,m,n);
end

%initialization初始化
[ii,jj] = find(H);             % subscript index to nonzero elements of H 将非零元素标上
indx = sub2ind(size(H),ii,jj); % linear index to nonzero elements of H
q0 = H * spdiags(f0(:),0,n,n);
sq0 = full(q0(indx)); 
sff0 = sq0;

q1 = H * spdiags(f1(:),0,n,n); 
sq1 = full(q1(indx));
sff1 = sq1;

%iterations
k=0;success = 0;max_iter =20;while ((success == 0) & (k < max_iter)),k = k+1;
%for iter=1:max_iter   
   
   %horizontal step第一步
   sdq = sq0 - sq1; sdq(find(sdq==0)) = 1e-20; % if   f0 = f1 = .5若先验概率一样的话，令sdq为一个很小很小的数
   dq = sparse(ii,jj,sdq,m,n);%以ii,jj为横纵坐标，用概率差sdq构成m*n的稀疏矩阵
   Pdq_v = full(real(exp(sum(spfun('log',dq),2)))); % this is ugly but works :)
                                                    %将dq在对数域求和，即将dq相乘，构成数组 Pdq_v                                                
   Pdq = spdiags(Pdq_v(:),0,m,m) * H;%更新了rij，以坐标形式6*12的矩阵
   sPdq = full(Pdq(indx));%变成了一列
   sr0 = (1+sPdq./sdq)./2; sr0(find(abs(sr0) < 1e-20)) = 1e-20;%将每个元素的rij都归一化，更新rij。概率译码的第一步完成
   sr1 = (1-sPdq./sdq)./2; sr1(find(abs(sr1) < 1e-20)) = 1e-20;%若比值非常小，则设为1e-20
   r0 = sparse(ii,jj,sr0,m,n);%将更新完的rij存入对应的位置，方便后面要用
   r1 = sparse(ii,jj,sr1,m,n);
   
   %vertical step第二步目的更新qij
   Pr0_v = full(real(exp(sum(spfun('log',r0),1))));
   Pr0 = H * spdiags(Pr0_v(:),0,n,n);
   sPr0 = full(Pr0(indx));
   Q0 = full(sum(sparse(ii,jj,sPr0.*sff0,m,n),1))';
   sq0 = sPr0.*sff0./sr0;
   
   Pr1_v = full(real(exp(sum(spfun('log',r1),1))));
   Pr1 = H * spdiags(Pr1_v(:),0,n,n);
   sPr1 = full(Pr1(indx)); 
   Q1 = full(sum(sparse(ii,jj,sPr1.*sff1,m,n),1))';
   sq1 = sPr1.*sff1./sr1;
   
   sqq = sq0+sq1;
   sq0 = sq0./sqq;
   sq1 = sq1./sqq;
   
   %tentative decoding伪后验概率（判决信息）
   QQ = Q0+Q1;
   Q0 = Q0./QQ;
   Q1 = Q1./QQ;
   
   x_hat = (sign(Q1-Q0)+1)/2;%判决输出
    if rem(H*x_hat,2) == 0, success = 1; 
        break;
    end
end
