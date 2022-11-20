% Def constants
K = 30;
r = 0.05;
% Necessary constraint
sigma = sqrt(2*r); 

% We set the time steps.
dtau = 0.001;
ds = 0.2;

tau_max = 5;

% Creates a vector of time units and x units.
tau = 0:dtau:tau_max;
s = 0:ds:5;
x = exp(s);
t = tau_max - tau;

% Counts the number of values to calculate in time(t) and space(x).
N = numel(tau);
M = numel(s);

% Creates an NxM matrix of zeros. 
u = zeros(N, M);
v_bar = zeros(N,M);

% Initial/Boundary conditions.
f = @(s) max(0, exp(s)-K);
%g = @(s) s*exp(t*r) - K;

u(1, :) = f(s);
u(N, end) = exp(r*tau(N) + s(end)) - K;
%u(end, :) = 0;
%u(:, end) = g(s);

% Applying finite difference approximation
for i = 1:N-1
    u(i, end) = exp(r*tau(i) + s(end)) - K;
    for j = 2:M-1
        u(end, j) = 0;
        u(i+1, j) = u(i, j) + (1/2)*(sigma^2)*(dtau/ds^2) * (u(i, j+1) -2*u(i, j) + u(i, j-1));
    end
end

for i = 1:N
    for j = 1:M
        v_bar(i,j) = u(i,j) * exp(-r*tau(i));
    end
end

% Graphs
waterfall(x,t,v_bar);
xlabel('Spot Price', 'FontSize', 15);
ylabel('Time', 'FontSize', 15);
zlabel('Option Value', 'FontSize', 15);
title('Black Scholes FD with K=30', 'FontSize', 20)
set(gca,'FontSize',20)