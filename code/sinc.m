function y = sinc(x)
%SINC - Compute the sinc function
%
%    y = SINC(x) is equal to sin(x)/x for nonzero x, and 1 for x = 0

% one at x = 0
y = ones(size(x));

% indices of non-zero x
iz = (x~=0);
x = x(iz);

y(iz) = sin(x) ./ x;
