function f = polar_design(biterrd,d)
% function f = polar_design(biterrd,d)
%   Design a polar code based on error rates of effective channels

% Sort into increasing order and compute cumulative sum
[SE,order] = sort(biterrd);
CSE = cumsum(SE);

% Find best frozen bits
k = sum(double(CSE<d));
f = zeros(1,length(biterrd));
f(order(1:k)) = 1/2;

