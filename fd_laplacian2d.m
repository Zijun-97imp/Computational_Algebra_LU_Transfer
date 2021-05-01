function A = fd_laplacian2d(h)

% A = FD_LAPLACIAN2D(h) constructs the finite difference matrix for the
% two-dimensional laplacian with Dirichlet boundary conditions on the unit
% square, meshed with step h in both directions. 
% Use sparse matrices the proper way.

n = 1/h - 1;                % Size of blocks in A
nnzT = 3*n-2;               % # nonzero entries in T
nnzA = n*nnzT + 2*n*(n-1);  % # nonzero entries of A
I = zeros(nnzA,1);  % container for row indices
J = zeros(nnzA,1);  % container for column indices
V = zeros(nnzA,1);  % container for nonzero values in A

% Nonzero entries due to the identity matrices on the off-diagonal blocks
s = 2*n*(n-1);
I(1:s) = [1:n*(n-1), n+1:n^2];
J(1:s) = [n+1:n^2, 1:n*(n-1)];
V(1:s) = -1;

% Diagonal entries
I(s+1:s+n^2) = 1:n^2;
J(s+1:s+n^2) = 1:n^2;
V(s+1:s+n^2) = 4;

% Entries above and below the main diagonal
s = s+n^2;
m = 2*(n-1);
r = 1:(n-1);
c = 2:n;
for i = 1:n
    I(s+1:s+m) = (i-1)*n + [r, c];
    J(s+1:s+m) = (i-1)*n + [c, r];
    V(s+1:s+m) = -1;
    s = s + m;
end

% Build sparse A
A = sparse(I,J,V,n^2,n^2);

