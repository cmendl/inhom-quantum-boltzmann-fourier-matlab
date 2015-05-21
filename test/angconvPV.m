function y = angconvPV(c1,c2,gamma)
%ANGCONVPV - trigonometric convolution using expansion by Bessel functions
%
%    y = ANGCONVPV(c1,c2,gamma)
%
%    Calculate the convolution
%    \int_0^{2*pi} d phi \int_0^{2*pi} d psi 1/cos(phi - psi) f(c1*cos(phi-alpha)) f(c2*cos(psi-beta))
%    with f(x) = sin(x) sinc(x), using expansion by Bessel J functions.
%    The integral only depends on gamma = alpha - beta.

if c1*c2 == 0
	y = 0;
	return;
end

besseljsum = @(j,z) (1-besselj(0,z) - 2*sum(besselj(2*(1:j),z)))/z;

% infinite sum cutoff
jmax = 50;

y = 0;
for j=0:jmax
	y = y + (-1)^j*cos((2*j+1)*gamma) * besseljsum(j,2*c1) * besseljsum(j,2*c2);
end
% scale factor
y = 2*(2*pi)^2 * y;
