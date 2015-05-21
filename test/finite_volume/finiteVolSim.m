function finiteVolSim()
%FINITEVOLSIM - Finite volume simulation example for u_t + A u_x = 0 with periodic boundary conditions
%
%	FINITEVOLSIM()

% eigenvalues of A
lambda = [-1,0.4,1.2];

% pseudo-random eigenvectors
R = [ 0.4133   -0.2381   -0.3634;
      0.2962   -0.1646    0.2212;
     -0.4013    0.1797   -0.3932];

% R*diag(lambda)*R^{-1}
A = R*diag(lambda)/R;

% spatial mesh width
h = 0.08;
% time step
k = 0.02;

% Courant number, Eq. (10.56), must be <= 1
courant = max(abs(eig(A)))*k/h;
fprintf('Courant number: %g (should be <= 1)\n',courant);

% spatial grid
xmax = 16;
x = 0:h:(xmax-h);
% time grid
tmax = 4;
t = 0:k:tmax;

% initial finite volume state
Uinit = zeros(length(A),length(x));
Uinit(1,:) = 3*max(1-abs(0.5*(x-1)),0);
Uinit(2,:) = 0.25*sin(pi/2*x);
Uinit(3,length(x)/2:end) = 1;	% jump
%Uinit(:, 1:end/2) = 1;

Ugod = cell(1,length(t));
Uslm = cell(1,length(t));
Ugod{1} = Uinit;
Uslm{1} = Uinit;
% reference version, represent in eigenbasis of A
Uref = cell(1,length(t));
for j=1:length(x)
	Uref{1}(:,j) = R\Uslm{1}(:,j);
end
% run simulation
for n=2:length(t)
	Ugod{n} =  godunovStep(h,k,A,Ugod{n-1});
	Uslm{n} = slopeLimStep(h,k,A,Uslm{n-1});
	% uncoupled 1D systems
	for p=1:length(A)
		Uref{n}(p,:) = slopeLimStep(h,k,lambda(p),Uref{n-1}(p,:));
	end
end
% transform back
for n=1:length(t)
	for j=1:length(x)
		Uref{n}(:,j) = R*Uref{n}(:,j);
	end
end

% visualize simulation
figure();
for n=1:length(t)
	plot(x,Ugod{n}(1,:),'-.b',x,Ugod{n}(2,:),'-.g',x,Ugod{n}(3,:),'-.r',...
		 x,Uslm{n}(1,:),'b',  x,Uslm{n}(2,:),'g',  x,Uslm{n}(3,:),'r',...
		 x,Uref{n}(1,:),'-.k',x,Uref{n}(2,:),'-.k',x,Uref{n}(3,:),'-.k');
	xlabel('x');
	ylabel('U');
	title(sprintf('dashed: Godunov, continuous: slope limiter, black dashed: reference\nt = %g',t(n)));
	xlim([0,xmax]);
	ylim([-15,15]);
	pause(0.1);
end
