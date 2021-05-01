% diff_2d.m
% ----------------
% Solve the 2D poisson equation
%
% -(T_xx + T_yy) = f
%
% with Dirichlet boundary conditions on the unit square using finite
% differences. The resulting linear system is solved using a direct method.

% Clean up
clear
close all

% ----------------------------------------------------------------------- %
% Setup
% ----------------------------------------------------------------------- %
% Number of gridpoints in each direction
n = 2^8;

% Function f
% i = 2;  j = 3;
% f = @(x,y) cos(i*pi*x.*y) - sin(j*pi*x.*y);
f = @(x,y) (100*cos(4*x.*y)+x).*(x-x.^2).*(y-y.^2);

% Mesh. Get y before x so order of nodes is as in the lecture notes
% after vectorization. Also initialize matrix for solution
h = 1/(n+1);
[y,x] = meshgrid(0:h:1,0:h:1);
Tmat = zeros(n+2,n+2);


% ----------------------------------------------------------------------- %
% Build the equation
% ----------------------------------------------------------------------- %
% Build system A*T = h^2*f. Evaluate f only at the internal points
A = fd_laplacian2d(h);
b = f(x(2:end-1,2:end-1),y(2:end-1,2:end-1));
b = h^2.*b(:);


% ----------------------------------------------------------------------- %
% Solve
% ----------------------------------------------------------------------- %
% % Solve using a sparse Cholesky factorization (matlab) and reshape into a
% % matrix for plotting. Also plot A and its Cholesky factor R
% tstart = tic;
% R = chol(A);
% chol_time = toc(tstart);
% Tvec = R\(R.'\b);
% runtime = toc(tstart);
% Tmat(2:end-1,2:end-1) = reshape(Tvec,n,n);

% Solve using a permuted sparse Cholesky factorization (matlab) and reshape
% into a matrix for plotting. Also plot the permuted A and its cholesky
% factor R.
tstart = tic;
[R,~,P] = chol(A);
chol_time = toc(tstart);
Tvec = P*( R\(R.'\(P.'*b)) );
runtime = toc(tstart);
Tmat(2:end-1,2:end-1) = reshape(Tvec,n,n);

% Display stats
fprintf('\n           n: %6i\n',n)
fprintf('     size(A): %6i\n',n^2)
fprintf('      nnz(A): %6.4e\n',nnz(A))
fprintf('      nnz(R): %6.4e\n',nnz(R))
fprintf('  Chol. time: %6.4f seconds\n',chol_time)
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
figure('WindowStyle','docked');
[~,h] = contourf(x,y,Tmat);
h.LevelStep = h.LevelStep/5;
h.LineStyle = 'none';
hold on; axis square; colorbar; colormap jet; caxis([-0.01,0.2])
title(sprintf('T(x,y), n = %i',n)); xlabel('x'); ylabel('y')