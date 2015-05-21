function err = I2int_test(N,J)
%I2INT_TEST - Compare I2 integral with reference implementation
%
%    err = I2INT_TEST()
%    err = I2INT_TEST(N,J)

% default parameters
if (nargin < 2)
	N = 8;
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
	psiR1{j} = crand(N);
	psiR2{j} = crand(N);
end

I2    = I2integral(f,g,h,psiR1,psiR2,weight);
I2ref = I2int_ref (f,g,h,psiR1,psiR2,weight);

% compare
err = norm(I2 - I2ref);
