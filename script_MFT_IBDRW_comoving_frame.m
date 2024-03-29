%--------------------------------------------------------------------------
% Inflationary theory of branching morphogenesis in the mouse salivary gland
% Bordeu I, Chatzeli L, and Simons BD (2022).
%--------------------------------------------------------------------------
% Numerical integration of the rescaled 2d mean-field theory of the 
% inflationary branching-delayed random walk (IBDRW) model, Eqs. (S15).
%
% See supplementary Note for details.
%
% The simulation shows a 1-d cut for a 2-dimensional simulation in the
% positive half-plane. The boundary conditions at x=0 and x=L are of 
% no-flux type, so that dA/dx=0 at the boundaries, and periodic elsewhere.
%
% For questions: ib443 (at) cam.ac.uk
%--------------------------------------------------------------------------
clear all;close all;clc;
%% Parameters -------------------------------------------------------------
% Mean-field parameters: 
rbranch = 1; % branching-to-elongation rate ratio
rexp = rbranch/8; % expansion-to-elongation rate ratio 
% integration parameters:
L = [201,1]; % system size
dt = 0.005; % time-step
dx = 0.1; % grid-spacing
N_iters = 15001; % Number of iterations
plot_times = 1:100:15001; % times at which to show density profiles.
% plot_times = [2000,4000,9000]; % used in the paper
%% Main -------------------------------------------------------------------
% Compute homogeneous steady-state solutions:
dim = 2;
h = dim*rexp; % 
a_plus = (dim-1)*(rbranch-h)*(h+rexp)/(rbranch*(rbranch*dim+dim-rbranch));
s_plus = (dim-1)*(rbranch-h)^2*(h+rexp)/(rexp*dim*rbranch*(rbranch*dim+dim-rbranch));
i_plus = (rbranch-h)*(h+rexp)/(rexp*rbranch*(rbranch*dim+dim-rbranch));
% Initialize system
Lx = L(1); Ly = L(2);
a = zeros(Ly,Lx); 
a(:,1) = a_plus/10; % initialise with a perturbation at A(x=0)
s = zeros(Ly,Lx);
i = zeros(Ly,Lx);
%
% Initialize figure
f = figure;
f.Position = [100 100 1700 600];
set(gcf,'color','w');
%
counter = 1; 
pos_front = [];
% run for N_iters timesteps
for niter = 1:N_iters
    % compute phi at every timestep
    phi = exp(rexp*dt*niter);
    % compute nabla^2 term for every field with appropriate boundary
    % conditions
    nabla2A = diffusion2d_9pt(a,dx);     
    % First order integration
    a = a + dt*(nabla2A/phi^2 + rbranch*a.*(1-(a + i + s)) + rexp*s - h*a);
    s = s + dt*(rbranch*a.*(a + i + s) - (h + rexp).*s);
    i = i + dt*(a - rexp*i*(dim-1));
    %
    % plot profiles, scaled by the homogeneous steady state solution
    if counter<=length(plot_times) && niter > plot_times(counter)
        counter = counter + 1;
        % Create figure object
        plot(a((Ly+1)/2,:)./a_plus,'k-');
        hold on
        plot(s((Ly+1)/2,:)./s_plus,'k-.');
        plot(i((Ly+1)/2,:)./i_plus,'k--');
        hold off
        refline(0,1)
        axis([0 Lx 0 3])
        xlabel('x'); ylabel('Normalized density')
        legend({'a: Active tips','s: Inactive tips',' i: Ducts','Homogeneous steady state (analytic solution)'},'Location','NorthEast'); legend boxoff
        pause(0.0001) 
    end
end

% FUNCTIONS ---------------------------------------------------------------
%
function nabla2A = diffusion2d_9pt(A,dx)
% input:
% A: NxN array (field)
% dx: spacing

% compute displacesd fields
Ap1x = circshift(A,-1,2);
Am1x = circshift(A,1,2);
Ap1y = circshift(A,-1,1);
Am1y = circshift(A,1,1);
Ap1p1 = circshift(Ap1x,-1,1);
Am1m1 = circshift(Am1x,1,1);
Ap1m1 = circshift(Ap1y,1,2);
Am1p1 = circshift(Am1y,-1,2);

% apply boundary conditions 
        Ap1x(:,end) = A(:,end);
        Am1x(:,1) = A(:,1);
        Ap1y(end,:) = A(end,:);
        Am1y(1,:) = A(1,:);
        
        Ap1p1(:,end) = A(:,end); Ap1p1(end,:) = A(end,:); 
        Am1m1(:,1) = A(:,1); Am1m1(1,:) = A(1,:);
        Ap1m1(end,:) = A(end,:); Ap1m1(:,1) = A(:,1);
        Am1p1(1,:) = A(1,:); Am1p1(:,end) = A(:,end);
% compute Laplacian
nabla2A = (Am1m1+ Ap1x + Am1x + Ap1p1 - 8.*A + Ap1m1 + Am1y + Ap1y + Am1p1)/(3*dx^2);
end % end diffusion2d_9pt
