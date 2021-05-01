function A = lu_direct_overwrite(A)

% A = LU_DIRECT_OVERWRITE(A) returns the LU factors of A, overwritten onto A.
%                  Precisely, the strictly lower triangular part of L is
%                  overwritten on the strictly lower triangular part of A,
%                  while the upper triangular factor U is overwritten on
%                  the upper triangular part of A.
%
% NOTE: This is an implementation of Algorithm 4.8 in the notes.

% Check inputs
[m,n] = size(A);
assert(m==n, 'Incorrect input size: A must be n-by-n')

% Compute LU factorization, overwrite factors on A
for k = 1:n-1
    if A(k,k) ~= 0
        % Use an "outer product" vectorized version of the loops over i and
        % j in algorithm 4.8 of the lecture notes. This is fast in matlab.
        A(k+1:n,k) = A(k+1:n,k)./A(k,k);
        A(k+1:n,k+1:n) = A(k+1:n,k+1:n) - A(k+1:n,k).*A(k,k+1:n);
    else
        error('Gaussian elimination failed: cannot divide by zero!')
    end
end