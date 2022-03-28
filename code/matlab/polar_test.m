% Setup code parameters
n = 12; N = 2^n;
e = 0.5; p = 0.10;
d = 0.1;
bec = 0;

% Compute the quality of all effective channels
if (bec)
  [biterrd] = polar_bec(n,e);
else
  [biterrd] = polar_bsc(n,p,1000);
end

% Design polar code
f = polar_design(biterrd,d);
A = (f==1/2);
k = sum(A);
rate = k/N

% Run a few sims to compare with union bound
M = 100;
biterr = zeros(1,M);
for i=1:M
  % Set frozen bits, add random data, and encode
  u = f;
  u(A) = rand(1,k)<0.5;
  x = polar_transform(u);

  % Transmit
  if (bec)
    y = x;
    y(rand(1,N)<e)=1/2;
  else
    y = zeros(1,N)+p;
    y(x==1) = 1-y(x==1);
    err = rand(1,N)<p;
    y(err) = 1-y(err);
  end

  % Decode and compute error rate for info bits
  [uhat,xhat] = polar_decode(y,f);
  biterr(i) = mean(uhat(A)~=u(A));
end

% Display average bit and block error rate
mean(biterr)
mean(biterr>0)
