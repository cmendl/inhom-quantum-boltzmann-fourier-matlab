function S = entropy(p)
%ENTROPY - Binary entropy function
%
%    S = ENTROPY(p)

% zero for p outside (0,1)
S = zeros(size(p));

% indices of valid p
iv = (p>0 & p<1);
p = p(iv);

S(iv) = -(p.*log(p) + (1-p).*log(1-p));
