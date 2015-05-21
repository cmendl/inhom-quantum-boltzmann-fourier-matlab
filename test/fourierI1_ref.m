function G = fourierI1_ref(N,L,R)
%FOURIERI1_REF - Reference implementation of the quadrature formula for the I1 integral for testing
%
%    G = FOURIERI1_REF(N,L,R)

kgrid1 = [0:N/2-1,-N/2:-1];
[kx1,ky1] = meshgrid(kgrid1,kgrid1);
% extended grid
kgrid2 = [0:N-1,-N:-1];
[kx2,ky2] = meshgrid(kgrid2,kgrid2);

G = zeros(N^2,N);
for l=1:numel(kx2)
	for m=1:numel(kx1)
		u = [kx2(l),ky2(l)];
		v = [kx1(m),ky1(m)];
		% polar angle difference
		dphi = atan2(u(2),u(1)) - atan2(v(2),v(1));
		% minus sign due to i^2
		G(l,m) = -R^2 * angconvPV(norm(u)*pi*R/(2*L),norm(v)*pi*R/(2*L),dphi);
	end
end
