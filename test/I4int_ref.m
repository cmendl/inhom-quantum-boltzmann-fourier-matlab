function I4 = I4int_ref(f,g,h,psiR1,psiR2,weight)
%I4INT_REF - Reference implementation of the I4 integral for testing
%
%    I4 = I4INT_REF(f,g,h,psiR1,psiR2,weight)
%
%    I4(v) = \int_{(R^2)^2} du du' \delta(u \cdot u') f(v + u) g(v + u + u') h(v)
%
%    f, g, h are N x N matrices storing the Fourier coefficients

% using the singular value decomposition
% H(\xi,\chi) = \sum_j weight(j) psiR1{j} psiR2{j}

N = size(f,1);
% N must be even
assert(mod(N,2) == 0);

toInd  = @(j) mod(j+N,N) + 1;
toInd2 = @(j) mod(j+2*N,2*N) + 1;

I4 = zeros(size(f));

for j=1:length(weight)

	% pointwise multiplication of two matrices
	g_psi = g .* psiR2{j};

	% \xi index
	for k1=-N/2:N/2-1
		for k2=-N/2:N/2-1

			% "\zeta loop"
			for l1=-N/2:N/2-1
				for l2=-N/2:N/2-1

					% "\eta loop"
					for m1=-N/2:N/2-1
						for m2=-N/2:N/2-1

							if (((k1-l1-m1) >= -N/2) && ((k1-l1-m1) < N/2))
								if (((k2-l2-m2) >= -N/2) && ((k2-l2-m2) < N/2))

									I4(toInd(k1),toInd(k2)) = I4(toInd(k1),toInd(k2)) + ...
										weight(j) * f(toInd(k1-l1-m1), ...
											  toInd(k2-l2-m2)) * g_psi(toInd(m1),toInd(m2)) * psiR1{j}(toInd2(k1-l1),toInd2(k2-l2)) * h(toInd(l1),toInd(l2));
								end
							end
						end
					end
				end
			end
		end
	end
end
