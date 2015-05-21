function I3 = I3int_ref(f,g,h,psiR1,psiR2,weight)
%I3INT_REF - Reference implementation of the I3 integral for testing
%
%    I3 = I3INT_REF(f,g,h,psiR1,psiR2,weight)
%
%    I3(v) = \int_{(R^2)^2} du du' \delta(u \cdot u') f(v + u) g(v + u') h(v + u + u')
%
%    f, g, h are N x N matrices storing the Fourier coefficients

% using the singular value decomposition
% H(\xi,\chi) = \sum_j weight(j) psiR1{j} psiR2{j}

N = size(f,1);
% N must be even
assert(mod(N,2) == 0);

toInd  = @(j) mod(j+N,N) + 1;
toInd2 = @(j) mod(j+2*N,2*N) + 1;

I3 = zeros(size(f));

for j=1:length(weight)
	% pointwise multiplication of two matrices
	f_psi = f .* psiR1{j};
	h_psi = h .* psiR1{j};

	t = zeros(2*size(f));

	% \eta index
	for k1=-N/2:N/2-1
		for k2=-N/2:N/2-1

			% "\zeta loop"
			for l1=-N/2:N/2-1
				for l2=-N/2:N/2-1

					m1 = toInd2(k1+l1);
					m2 = toInd2(k2+l2);

					t(m1,m2) = t(m1,m2) + g(toInd(k1),toInd(k2)) ...
						* psiR2{j}(m1,m2) * h_psi(toInd(l1),toInd(l2));
				end
			end
		end
	end

	% \xi index
	for k1=-N/2:N/2-1
		for k2=-N/2:N/2-1

			% "\chi loop"
			for l1=-N/2:N/2-1
				for l2=-N/2:N/2-1

					I3(toInd(k1),toInd(k2)) = I3(toInd(k1),toInd(k2)) +...
						weight(j) * t(toInd2(k1-l1),toInd2(k2-l2)) * f_psi(toInd(l1),toInd(l2));
				end
			end
		end
	end
end
