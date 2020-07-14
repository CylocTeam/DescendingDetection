function [y] = TwoPoleLPF(t, x, fres, damping )
if isempty(t)
    t = 1:size(x,1);
end
if fres == 0
    y=x;  % no bandwidth
    return
end
% w0 = 1 / tau;
w0 = 2* pi * fres;
R = damping * w0;
wd = w0 * sqrt(1 - damping ^ 2);
N = length(t);

FLIPPED = false;
if length(size(x))== 1
    x = x(:); 
end                 % make this a row vector
[l, D] = size(x);
    if D == N
        x = x';  % want x to be long along the zero dimension
        D = l;  % D is the dimension of x. I want x to be Nx3 for example
        FLIPPED = True;
    end
  
y = zeros(size(x));
y(1, :) = x(1, :);
y0 = x(1, :);
yDot0 = zeros(1, D);
for n = 2 : N
    dt = t(n) - t(n-1);
    x0 = x(n, :);
    A2 = y0 - x0;
    A1 = (yDot0 + R * A2) / wd;
    sinTerm = sin(wd * dt);
    cosTerm = cos(wd * dt);
    expTerm = exp(-R * dt);
    y0 = x0 + expTerm * (A1 * sinTerm + A2 * cosTerm);
    yDot0 = expTerm * ((wd * A1 - R * A2) * cosTerm - (wd * A2 + R * A1) * sinTerm);
    y(n, :) = y0;
end
if FLIPPED
    y = y';
    return
end
end
