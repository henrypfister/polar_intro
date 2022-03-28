function [u,x] = polar_decode(y,f)
%function [u,x] = polar_decode(y,f)
% y = channel observations in output order
% x = output hard decisions in output order
% f = input a priori probabilities in input order
% u = input hard decisions in input order

  % Recurse down to length 1
  N = length(y);
  if (N==1)
    % If information bit (i.e., f=1/2 for P1 domain)
    if (f==1/2)
      % Make hard decision based on observation
      x = (1-sign(1-2*y))/2; u = x;
    else
      % Use frozen bit
      u = (1-sign(1-2*y+1e-10*(rand(1,N)-0.5)))/2; x = f;
      %x = f; u = x;
    end
  else
    % Compute soft mapping back one stage
    u1est = cnop(y(1:2:end),y(2:2:end));

    % R_N^T maps u1est to top polar code
    [uhat1,u1hardprev] = polar_decode(u1est,f(1:(N/2)));

    % Using u1est and x1hard, we can estimate u2
    u2est = vnop(cnop(u1hardprev,y(1:2:end)),y(2:2:end));

    % R_N^T maps u2est to bottom polar code
    [uhat2,u2hardprev] = polar_decode(u2est,f((N/2+1):end));

    % Pass u decisions up and interleave x1,x2 hard decisions
    u = [uhat1 uhat2];
    x1 = cnop(u1hardprev,u2hardprev);
    x2 = u2hardprev;
    x = reshape([x1; x2],1,[]);
  end
return

% Check-node operation in P1 domain
function z=cnop(w1,w2)
  z = w1.*(1-w2) + w2.*(1-w1);
return

% Bit-node operation in P1 domain
function z=vnop(w1,w2)
  z = w1.*w2 ./ (w1.*w2 + (1-w1).*(1-w2));
return

