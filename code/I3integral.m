function I3 = I3integral(f,g,h,psiR1,psiR2,weight)
%I3INTEGRAL - Evaluate the I3 integral
%
%    I3 = I3INTEGRAL(f,g,h,psiR1,psiR2,weight)
%
%    I3(v) = \int_{(R^2)^2} du du' \delta(u \cdot u') f(v + u) g(v + u') h(v + u + u')
%
%    f, g, h are N x N matrices storing the Fourier coefficients with frequencies [-N/2,N/2-1]^2,
%    with the first entry corresponding to freqency zero, lineary increasing modulo N

I3 = zeros(2*size(f));

N = length(f);
IX = [1:N/2,3*N/2+1:2*N];

% pre-compute FFT of g2
g2 = zeros(2*N);
g2(IX,IX) = g;
fft_g2 = fft2(g2);

for j=1:length(weight)
	% pointwise multiplication of two matrices
	f_psi = f .* psiR1{j};
	h_psi = h .* psiR1{j};

	% zero padding for h_psi
	h_psi2 = zeros(2*N);
	h_psi2(IX,IX) = h_psi;

	% convolution of g2 and h_psi2
	t = ifft2(fft_g2.*fft2(h_psi2));

	% pointwise multiplication with psiR2{j}
	t = t .* psiR2{j};

	% zero padding for f_psi
	f_psi2 = zeros(2*N);
	f_psi2(IX,IX) = f_psi;

	% convolution with f_psi
	I3 = I3 + weight(j) * ifft2(fft2(f_psi2).*fft2(t));
end

% extract [-N/2:N/2-1] part
I3 = I3(IX,IX);
