function B = fd_advection2d(U,V)

% A = fd_advection2d(h) constructs the finite difference matrix for the
% two-dimensional advection term with Dirichlet boundary conditions on the
% unit square, meshed with step h in both directions.

% number of gridpoints in interior of the domain (T=0 on boundaries!)
n = size(U,1)-2;

% Initialize variables with correct size
nnzB = (n-2)^2*4 + 4*2 + 4*(n-2)*3;
I = zeros(nnzB,1);
J = zeros(nnzB,1);
VALS = zeros(nnzB,1);

% Loop over all gridpoints
shift = 0;
for k = 1:n^2
    
    % Find i,j indices
    j = fix((k-1)./n) + 1;
    i = k - n*(j-1);
    
    % What case do we have?
    if i==1 && j==1
        % bottom left corner
        J(shift+1:shift+2) = [n*(j-1)+(i+1); n*j+i];
        I(shift+1:shift+2) = [k; k];
        VALS(shift+1:shift+2) = [U(i+1,j); V(i,j+1)];
        shift = shift + 2;
        
    elseif i==n && j==1
        % bottom right corner
        J(shift+1:shift+2) = [i-1; n*j+i];
        I(shift+1:shift+2) = [k; k];
        VALS(shift+1:shift+2) = [-U(i-1,j); V(i,j+1)];
        shift = shift + 2;
        
    elseif i==1 && j==n
        % top left corner
        J(shift+1:shift+2) = [n*(j-1)+(i+1); n*(j-2)+i];
        I(shift+1:shift+2) = [k; k];
        VALS(shift+1:shift+2) = [U(i+1,j); -V(i,j-1)];
        shift = shift + 2;
        
    elseif i==n && j==n
        % top right corner
        J(shift+1:shift+2) = [n*(j-1)+(i-1); n*(j-2)+i];
        I(shift+1:shift+2) = [k; k];
        VALS(shift+1:shift+2) = [-U(i-1,j); -V(i,j-1)];
        shift = shift + 2;
        
    elseif i==1 % and 1<j<n
        % Left boundary, no corners
        J(shift+1:shift+3) = [n*(j-1)+(i+1); n*j+i; n*(j-2)+i];
        I(shift+1:shift+3) = [k; k; k];
        VALS(shift+1:shift+3) = [U(i+1,j); V(i,j+1); -V(i,j-1)];
        shift = shift + 3;
        
    elseif i==n % and 1<j<n
        % right boundary, no corners
        J(shift+1:shift+3) = [n*(j-1)+(i-1); n*j+i; n*(j-2)+i];
        I(shift+1:shift+3) = [k; k; k];
        VALS(shift+1:shift+3) = [-U(i-1,j); V(i,j+1); -V(i,j-1)];
        shift = shift + 3;
        
    elseif j==1 % and 1<i<n
        % bottom boundary, no corners
        J(shift+1:shift+3) = [n*(j-1)+(i+1); n*(j-1)+(i-1); n*j+i];
        I(shift+1:shift+3) = [k; k; k];
        VALS(shift+1:shift+3) = [U(i+1,j); -U(i-1,j); V(i,j+1)];
        shift = shift + 3;
        
    elseif j==n % and 1<i<n
        % top boundary, no corners
        J(shift+1:shift+3) = [n*(j-1)+(i+1); n*(j-1)+(i-1); n*(j-2)+i];
        I(shift+1:shift+3) = [k; k; k];
        VALS(shift+1:shift+3) = [U(i+1,j); -U(i-1,j); -V(i,j-1)];
        shift = shift + 3;
        
    else
        % Internal point
        J(shift+1:shift+4) = [n*(j-1)+(i+1); n*(j-1)+(i-1); n*j+i; n*(j-2)+i];
        I(shift+1:shift+4) = [k; k; k; k];
        VALS(shift+1:shift+4) = [U(i+1,j); -U(i-1,j); V(i,j+1); -V(i,j-1)];
        shift = shift + 4;
        
    end
end

% Build sparse A
B = sparse(I,J,VALS./2,n^2,n^2);