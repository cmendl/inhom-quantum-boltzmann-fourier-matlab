function I4 = I4integral(f,g,h,psiR1,psiR2,weight)
%I4INTEGRAL - Evaluate the I4 (or I1) integral
%
%    I4 = I4INTEGRAL(f,g,h,psiR1,psiR2,weight)
%
%    I4(v) = \int_{(R^2)^2} du du' \delta(u \cdot u') f(v + u) g(v + u + u') h(v)
%
%    f, g, h are N x N matrices storing the Fourier coefficients with frequencies [-N/2,N/2-1]^2,
%    with the first entry corresponding to freqency zero, lineary increasing modulo N

% using the singular value decomposition
% H(\xi,\chi) = \sum_j weight(j) psiR1{j} psiR2{j}

% extended version
I4 = zeros(2*size(f));

N = length(f);
IX = [1:N/2,3*N/2+1:2*N];

% pre-compute FFT of f
f2 = zeros(2*N);
f2(IX,IX) = f;
fft_f2 = fft2(f2);

% zero padding for h and FFT
h2 = zeros(2*N);
h2(IX, IX) = h;
fft_h2 = fft2(h2);

for j=1:length(weight)

	% pointwise multiplication of two matrices
	g_psi = g .* psiR2{j};
	g_psi2 = zeros(2*N);
	g_psi2(IX,IX) = g_psi;

	% convolution of f with g_psi
	fg2 = ifft2(fft_f2.*fft2(g_psi2));

	% pointwise multiplication with psiR1
	fg_psi = fg2 .* psiR1{j};

	% convolution with h
	I4 = I4 + weight(j) * ifft2(fft2(fg_psi).*fft_h2);
end

% extract [-N/2:N/2-1] part

I4 = I4(IX,IX);
