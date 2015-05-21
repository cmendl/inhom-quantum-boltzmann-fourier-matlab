function err = conv_test(N)
%CONV_TEST - Test file concerning convolution in 1D
%
%    err = CONV_TEST(N)

% N must be even
assert(mod(N,2) == 0);

crand = @(n) (rand(n)-0.5) + 1i*(rand(n)-0.5);

% random numbers
f = crand([1,N]);
g = crand([1,N]);

% The first entry in f and g stores Fourier coefficient for frequency zero,
% and negative frequencies are shifted by N in index space

toIndex = @(j) mod(j+N,N) + 1;

% reference implementation
h_ref = zeros(size(f));
for j=-N/2:N/2-1
	for k=-N/2:N/2-1
		h_ref(toIndex(j)) = h_ref(toIndex(j)) + f(toIndex(k))*g(toIndex(j-k));
	end
end

% implementation using FFT
h = ifft(fft(f) .* fft(g));

% compare
err = norm(h - h_ref);
