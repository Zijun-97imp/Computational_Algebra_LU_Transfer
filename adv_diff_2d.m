% adv_diff_2d.m
% ----------------
% Solve a 2D steady advection-diffusion equation
%
% -(T_xx + T_yy) + (uT)_x + (vT)_y = f
%
% with Dirichlet boundary conditions on the unit square using finite
% differences. The resulting linear system is solved using a direct method.

% Clean up
clear

% ----------------------------------------------------------------------- %
% Setup
% ----------------------------------------------------------------------- %
% Number of gridpoints in each direction
n = 2^8;

% source term
% i = 2;  j = 3;
% f = @(x,y) cos(i*pi*x.*y) - sin(j*pi*x.*y);
f = @(x,y) (100*cos(4*x.*y)+x).*(x-x.^2).*(y-y.^2);

% Advection velocities in the two directions: rolls with free-slip boundary
% conditions
l = 3; m = 3; Umag = 10;
u = @(x,y) Umag*m*sin(l*pi*x).*cos(m*pi*y);
v = @(x,y) -(Umag*l).*cos(l*pi*x).*sin(m*pi*y);

% Mesh. Get y before x so order of nodes is as in the lecture notes
% after vectorization. Also initialize matrix to store the solution.
% Finally, calculate the flow velocities U and V
h = 1/(n+1);
[y,x] = meshgrid(0:h:1,0:h:1);
Tmat = zeros(n+2,n+2);
Uvel = u(x,y);
Vvel = v(x,y);


% ----------------------------------------------------------------------- %
% Build the equation
% ----------------------------------------------------------------------- %
% Build the *negative* Laplacian matrix (multiplied by h^2)
DIFF = fd_laplacian2d(h);

% Build the advection matrix (multiplied by h). Need to multiply it by h,
% which gives advection matrix * h^2 in the discretized advection-diffusion
% equation (the h^2 term comes from the laplacian)
ADV = fd_advection2d(Uvel,Vvel);
ADV = h.*ADV;

% Make coefficient matrix M: diffusion + advection. This is the left-hand
% side of the equation. It is NOT symmetric!
A = DIFF + ADV;

% Right-hand side: evaluate f (only at the internal points) and make it a
% vector, ready for solution.
b = f(x(2:end-1,2:end-1),y(2:end-1,2:end-1));
b = h^2.*b(:);


% ----------------------------------------------------------------------- %
% Solve
% ----------------------------------------------------------------------- %
% Solve using a sparse LU factorization with partial pivoting (matlab
% built-in) and reshape into a matrix for plotting.
tstart = tic;
[L,U,p] = lu(A,'vector');
lu_time = toc(tstart);
Tvec = U\(L\b(p));
runtime = toc(tstart);
Tmat(2:end-1,2:end-1) = reshape(Tvec,n,n);

% Display stats
fprintf('\n           n: %6i\n',n)
fprintf('     size(A): %6i\n',n^2)
fprintf('      nnz(A): %6.4e\n',nnz(A))
fprintf('      nnz(L): %6.4e\n',nnz(L))
fprintf('      nnz(U): %6.4e\n',nnz(U))
fprintf('     LU time: %6.4f seconds\n',lu_time)
fprintf('  Total time: %6.4f seconds\n\n',runtime)

% ----------------------------------------------------------------------- %
% Plot results
% ----------------------------------------------------------------------- %
% Plot f
figure('WindowStyle','docked')
[~,h] = contourf(x,y,f(x,y));
h.LevelStep = h.LevelStep/3;
h.LineStyle = 'none';
axis square; colorbar; colormap jet
title('f(x,y)'); xlabel('x'); ylabel('y')

% Plot solution
figure('WindowStyle','docked')
[~,h] = contourf(x,y,Tmat);
h.LevelStep = h.LevelStep/5;
h.LineStyle = 'none';
hold on

% Overlay coarse flow field: quiver plot
[yrough, xrough] = meshgrid(0:0.05:1,0:0.05:1);
quiver(xrough,yrough,u(xrough,yrough),v(xrough,yrough),'w','linewidth',1);
axis([0 1 0 1]); axis square; colorbar; colormap jet; caxis([-0.01,0.2])
title(sprintf('T(x,y), n = %i',n)); xlabel('x'); ylabel('y')