function A = random_A(n,bands)

% A = RANDA(n) generates a random n-by-n matrix A for which basic Gaussian
%              elimination does not fail.
%
% A = RANDA(n,[ml,mu]) generates a random n-by-n banded matrix A for which basic 
%              Gaussian elimination does not fail. A has upper bandwidth mu
%              and lower bandwidth ml.

if nargin < 2; bands = [n-1,n-1]; end

% Trick: build A from invertible LU factors
% L = n*rand(n); L = eye(n) + triu(tril(L,-1),-bands(1));
% U = n*rand(n); U = eye(n) + tril(triu(U),bands(2));
L = randi(2,n); L = eye(n) + triu(tril(L,-1),-bands(1));
U = randi(2,n); U = eye(n) + tril(triu(U),bands(2));
A = L*U;