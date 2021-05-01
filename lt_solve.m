function x = lt_solve(L,b)

% x = lt_solve(L,b) solves the lower triangular system L*x=b using forward
% substitution. Only the lower triangular part of L is used.

% Preliminary stuff: size of L and preallocate solution (column vector)
% for speed
nL = size(L,1);
x = zeros(nL,1);

% Compute (no error checks)
x(1) =  b(1)/L(1,1);
for i = 2:nL
    x(i) =  ( b(i) - L(i,1:i-1)*x(1:i-1) ) / L(i,i);
end