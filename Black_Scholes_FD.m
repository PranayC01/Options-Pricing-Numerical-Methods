% We set the time steps.
dt = 0.0001;
ds = 0.02;

% Creates a vector of time units and x units.
t = 0:dt:5;
s = 0:ds:5;

% Counts the number of values to calculate in time(t) and space(x).
N = numel(t);
M = numel(s);

% Creates an NxM matrix of zeros. 
u = zeros(N, M);

% Def constants
K = 1;
r = 0.01;
sigma = sqrt(2*r);

% Initial/Boundary conditions.
f = @(s) max(0, exp(s)-K);
%g = @(s) s*exp(t*r) - K;

u(1, :) = f(s);
%u(end, :) = 0;
%u(:, end) = g(s);

% Applying finite difference approximation
for i = 1:N-1
    for j = 2:M-1
        u(end, j) = 0;
        u(i+1, j) = u(i, j) + (1/2)*(sigma^2)*(dt/ds^2) * (u(i, j+1) -2*u(i, j) + u(i, j-1));
    end
end

% Graphs
waterfall(s,t,u);
xlabel('s', 'FontSize', 15);
ylabel('t', 'FontSize', 15);
zlabel('u', 'FontSize', 15);
title('Black Scholes FD', 'FontSize', 30)
set(gca,'FontSize',20)