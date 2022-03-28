function E = polar_bec(n,e)
% function E = polar_bec(n,e)
%   Compute effective-channel erasure rates for polar code of length N=2^n on BEC(e)
E = e;
for i=1:n
  % Interleave updates to keep in polar decoding order
  E = reshape([1-(1-E).*(1-E); E.*E],1,[]);
end
