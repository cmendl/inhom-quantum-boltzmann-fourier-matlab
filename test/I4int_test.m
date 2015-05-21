function err = I4int_test(N,J)
%I4INT_TEST - Compare I4 integral with reference implementation
%
%    err = I4INT_TEST()
%    err = I4INT_TEST(N,J)

% default parameters
if (nargin < 2)
	N = 6;
	J = 7;
end

% N must be even
assert(mod(N,2) == 0);

crand = @(n) (rand(n)-0.5) + 1i*(rand(n)-0.5);

% random matrix entries
f = crand(N);
g = crand(N);
h = crand(N);

weight = rand(1,J);
psiR1 = cell(1,J);
psiR2 = cell(1,J);
for j=1:J
	psiR1{j} = crand(2*N);
	psiR2{j} = crand(N);
end

I4    = I4integral(f,g,h,psiR1,psiR2,weight);
I4ref = I4int_ref (f,g,h,psiR1,psiR2,weight);

% compare
err = norm(I4 - I4ref);
