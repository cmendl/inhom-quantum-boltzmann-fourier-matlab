function I2 = I2int_ref(f,g,h,psiR1,psiR2,weight)
%I2INT_REF - Reference implementation of the I2 integral for testing
%
%    I2 = I2INT_REF(f,g,h,psiR1,psiR2,weight)
%
%    I2(v) = \int_{(R^2)^2} du du' \delta(u \cdot u') f(v + u) g(v + u') h(v)
%
%    f, g, h are N x N matrices storing the Fourier coefficients

% using the singular value decomposition
% H(\xi,\chi) = \sum_j weight(j) psiR1{j} psiR2{j}

N = size(f,1);
% N must be even
assert(mod(N,2) == 0);

toInd  = @(j) mod(j+N,N) + 1;
toInd2 = @(j) mod(j+2*N, 2*N) + 1;

I2 = zeros(size(f));

for j=1:length(weight)

	% pointwise multiplication of two matrices,
	% implements multiplication with H(\xi, \chi)
	f_psi = f .* psiR1{j};
	g_psi = g .* psiR2{j};

	% convolution of f_psi and g_psi
	t = zeros(2*size(f));

	% "\chi loop"
	for k1=-N/2:N/2-1
		for k2=-N/2:N/2-1
			% "\eta loop"
			for l1=-N/2:N/2-1
				for l2=-N/2:N/2-1
					% index for t
					m1 = toInd2(k1+l1);
					m2 = toInd2(k2+l2);
					t(m1,m2) = t(m1,m2) + f_psi(toInd(k1),toInd(k2)) * g_psi(toInd(l1),toInd(l2));
				end
			end
		end
	end

	% convolution with h

	% \xi index
	for k1=-N/2:N/2-1
		for k2=-N/2:N/2-1
			% "\zeta loop"
			for l1=-N/2:N/2-1
				for l2=-N/2:N/2-1
					I2(toInd(k1),toInd(k2)) = I2(toInd(k1),toInd(k2)) + ...
						weight(j) * t(toInd2(k1-l1),toInd2(k2-l2)) * h(toInd(l1),toInd(l2));
				end
			end
		end
	end
end
