% 2x2 example of splitting methods
clear
close all

% Define system & convergence tolerance
A = [3, -1; -1, 3];
b = [2; -1];
TOL = 1e-12;

% Define the four cases (M only, find N=A-M)
M(:,:,1) = [3, 0; 0, 3];
M(:,:,2) = [3, 0; -1, 3];
M(:,:,3) = [0, -1; -1, 0];
M(:,:,4) = [0, -1; -1, 3];

% Loop over different test cases. Stop and wait for user to continue before
% continuing to next case
for i = 1:4
    
    % Calculate stuff
    N = A - M(:,:,i);               % Get N
    x = [0; 0];                     % Initial guess
    res = norm(A*x-b,Inf);          % Initial residual
    k = 0;                          % iteration number
    
    % Loop - the iterative method
    % Keep log of the residual error
    while res(k+1) > TOL && k < 10
        k = k + 1;
        x = M(:,:,i)\(b-N*x);
        res(k+1) = norm(A*x-b,Inf);
    end
    
    % Plot results in semilog scale
    semilogy(0:k,res,'o-','DisplayName',sprintf('Method %i',i),'LineWidth',1.25,'markersize',8);
    hold on
    if i==1; ll = legend('toggle'); end
    ll.Interpreter = 'latex';
    ll.Location = 'eastoutside';
    ll.FontSize = 16;
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    ax.FontSize = 16;
    ax.YTick = [1e-10 1e-5 1 1e5 1e10];
    xlabel('$k$','Interpreter','latex','FontSize',16)
    ylabel('$\|A\mathbf{x}_k - \mathbf{b}\|_\infty$','Interpreter','latex','FontSize',16)
    drawnow; pause();
end