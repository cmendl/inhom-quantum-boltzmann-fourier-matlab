function finiteVolSimMaxwell()
%FINITEVOLSIMMAXWELL - Finite volume simulation example for u_t + A u_x = 0 with Maxwell boundary conditions
%
%	FINITEVOLSIMMAXWELL()

% eigenvalues of A
lambda = [1,0.5,-0.5,-1];

% pseudo-random eigenvectors
R = [ 0.0466    0.1790    0.2093   -0.0499;
     -0.0743    0.1358   -0.2638   -0.0413;
      0.1444    0.4452   -0.3806    0.1619;
      0.1476   -0.2911    0.1073    0.2703];

% R*diag(lambda)*R^{-1}
A = R*diag(lambda)/R;

% Eq. (10.61)
Apos = R*diag(max(lambda,0))/R;
Aneg = R*diag(min(lambda,0))/R;

% spatial mesh width
h = 0.01;
% time step
dt = 0.005;

% accommodation coefficient for Maxwell reflection operator
alpha = 0.4;

% Courant number, Eq. (10.56), must be <= 1
courant = max(abs(eig(A)))*dt/h;
fprintf('Courant number: %g (should be <= 1)\n',courant);

% spatial grid
xmax = 1;
x = 0:h:xmax;
% time grid
tmax = 0.4;
t = 0:dt:tmax;

% initial finite volume state
Uinit = zeros(length(A),length(x));
Uinit(1,:) = x+0.4;
Uinit(2,:) = cos(pi*x);
Uinit(3,:) = -1;
Uinit(3,ceil(end/2):end) = 1;	% jump
Uinit(4,:) = erf(5*(x-0.5));

% incoming state from the left
UmaxwL = [0.2;-pi/4;0;0];		% non-zero for velocities > 0

% incoming state from the right
UmaxwR = [0;0;exp(-1);-0.8];	% non-zero for velocities < 0

% normalize prescribed boundary states
UmaxwL = UmaxwL / sum( Apos*UmaxwL);
UmaxwR = UmaxwR / sum(-Aneg*UmaxwR);

% run simulation
Ugod = cell(1,length(t));
Uslm = cell(1,length(t));
Ugod{1} = Uinit;
Uslm{1} = Uinit;
for n=2:length(t)
	Ugod{n} =  godunovStepMaxwell(h,dt,A,alpha,UmaxwL,UmaxwR,Ugod{n-1});

	% left and right flux
	fluxL = sum(-Aneg*Uslm{n-1}(:,1));		% minus sign from normal pointing in negative x direction
	fluxR = sum( Apos*Uslm{n-1}(:,end));
	Uslm{n} = slopeLimStepMaxwell(h,dt,A,alpha,fluxL,fluxR,UmaxwL,UmaxwR,Uslm{n-1});
end

% visualize simulation
for n=1:length(t)
	plot(x,Ugod{n}(1,:),'-.b',x,Ugod{n}(2,:),'-.g',x,Ugod{n}(3,:),'-.r',x,Ugod{n}(4,:),'-.m',...
		 x,Uslm{n}(1,:),'b',  x,Uslm{n}(2,:),'g',  x,Uslm{n}(3,:),'r',  x,Uslm{n}(4,:),'m');
	xlabel('x');
	ylabel('U');
	title(sprintf('dashed: Godunov, continuous: slope limiter\nt = %g',t(n)));
	xlim([0,xmax]);
	ylim([-10,10]);
	pause(0.1);
end
