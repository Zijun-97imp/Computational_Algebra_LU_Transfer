function x = ut_solve(U,b)

% x = ut_solve(U,b) solves the upper triangular system U*x=b using backward
% substitution. Only the upper triangular part of L is used.

% Preliminary stuff: size of U and preallocate solution (column vector)
% for speed
nU = size(U,1);
x = zeros(nU,1);

% Compute (no error checks)
x(nU) =  b(nU)/U(nU,nU);
for i = nU-1:-1:1
    x(i) =  ( b(i) - U(i,i+1:nU)*x(i+1:nU) ) / U(i,i);
end