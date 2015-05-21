function I2 = I2integral(f,g,h,psiR1,psiR2,weight)
%I2INTEGRAL - Evaluate the I2 integral
%
%    I2 = I2INTEGRAL(f,g,h,psiR1,psiR2,weight)
%
%    I2(v) = \int_{(R^2)^2} du du' \delta(u \cdot u') f(v + u) g(v + u') h(v)
%
%    f, g, h are N x N matrices storing the Fourier coefficients with frequencies [-N/2,N/2-1]^2,
%    with the first entry corresponding to freqency zero, lineary increasing modulo N

% using the singular value decomposition
% H(\xi,\chi) = \sum_j weight(j) psiR1{j} psiR2{j}

I2 = zeros(2*size(f));

N = length(f);
IX = [1:N/2,3*N/2+1:2*N];

% pre-compute FFT of h
h2 = zeros(2*N);
h2(IX,IX) = h;
fft_h2 = fft2(h2);

for j=1:length(weight)
	% pointwise multiplication of two matrices
	f_psi = f .* psiR1{j};
	g_psi = g .* psiR2{j};

	% extend
	f_psi2 = zeros(2*N);
	f_psi2(IX,IX) = f_psi;

	g_psi2 = zeros(2*N);
	g_psi2(IX,IX) = g_psi;

	% convolution of f_psi, g_psi and h
	I2 = I2 + weight(j) * ifft2(fft2(f_psi2).*fft2(g_psi2).*fft_h2);
end

% extract [-N/2:N/2-1] part
I2 = I2(IX,IX);
