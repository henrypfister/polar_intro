function [biterrd] = polar_bsc(n,p,M)
% function [biterrd] = polar_bsc(n,p,M)
%   Send M blocks for Monte Carlo estimate of length N=2^n polar code on BSC(p)

% Setup parameters
N = 2^n;
f = zeros(1,N);
biterrd = zeros(1,N);

% Monte Carlo evaluation of error probability
for i=1:M
  % Transmit all-zero codeword through BSC(p)
  y = zeros(1,N)+p;
  y(rand(1,N)<p)=1-p;
  % Decode received vector using all-zero frozen vector
  [uhat,xhat] = polar_decode(y,f);
  biterrd = biterrd + uhat;
end
biterrd = biterrd/M;

