function [waveform]=bpsk(bitseq)
%   maps binary 0 to -1 (and 1 to +1, the function actually does not do
%   this , as these values are already at (+)1  )


waveform=ones(1,length(bitseq));
waveform(find(bitseq==0))=-1;
