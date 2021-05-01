function [L, U] = lu_direct(A)

% [L, U] = LU_DIRECT(A) returns the LU factors of A, where L is a lower
%          triangular matrix with diagonal entries equal to 1 and U is the
%          upper triangular matrix obtained via Gaussian elimination.
%
% NOTE: This is an implementation of Algorithm 4.7 in the notes.

% Check inputs
[m,n] = size(A);
assert(m==n, 'Incorrect input size: A must be n-by-n')

% Compute
U = A;          % Initalize U
L = eye(n);     % Initialize L
for k = 1:n-1
    if U(k,k) ~= 0
        % Use an "outer product" vectorized version of the loops over i and
        % j in algorithm 4.7 of the lecture notes. This is fast in matlab.
        L(k+1:n,k) = U(k+1:n,k)/U(k,k);
        U(k+1:n,k:n) = U(k+1:n,k:n) - L(k+1:n,k).*U(k,k:n);
    else
        error('Gaussian elimination failed: cannot divide by zero!')
    end
end