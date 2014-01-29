% FORCES - Fast interior point code generation for multistage problems.
% Copyright (C) 2011-12 Alexander Domahidi [domahidi@control.ee.ethz.ch],
% Automatic Control Laboratory, ETH Zurich.
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


% MPC with increment input limits and soft outout bounds
% 
%  min   xN'*P*xN + sum_{i=1}^{N-1} xi'*Q*xi + ui'*R*ui + dui'*S*dui
% xi,ui
%       s.t. x1 = x, dui = ui - u(i-1)
%            x_i+1 = A*xi + B*ui  for i = 1...N-1
%            xmin-eli <= xi <= xmax+eui   for i = 1...N
%            umin <= ui <= umax   for i = 1...N
%             eli,eui >= 0 for i = 1...N
%
% and P is solution of Ricatti eqn. from LQR problem

% make sure you have FORCES on your path
addpath('FORCES');

%% system
nx = 4;
nu = 2;
A = [0.0211, 0,       0,      0;
     0.1062, 0.4266,  0,      0;
     0,      0,       0.2837, 0;
     0.1012, -0.6688, 0.2893, 0.4266];
B = [0.6462, 0.6462;
     0.2800, 0.2800;
     1.5237, -0.7391;
     0.9929, 0.1507];
C = eye(nx);
Q = eye(nx);
R = eye(nu);
S = eye(nu);

umin = [-2; -2]; umax = [2; 2];
xmin = [0; 0; 0; 0]; xmax = [10; 10; 10; 10];
dumax = [0.1; 0.1];

%% enhanced system
Qu = 0.1*eye(nx); Ql = 0.1*eye(nx);
Aa = [A, zeros(nx,nu); zeros(nu,nx), zeros(nu)];
Ba = [B, zeros(nx), zeros(nx); eye(nu), zeros(nu,nx), zeros(nu,nx)];
Q = [C'*Q*C, zeros(nx,nu); zeros(nu,nx), S];
R = [R+S, zeros(nu,nx), zeros(nu,nx); zeros(nx,nu), Qu, zeros(nx); zeros(nx,nu), zeros(nx), Ql];
M = [zeros(nx,nu),zeros(nx), zeros(nx); -S, zeros(nu,nx), zeros(nu,nx)];
umaxa = [umax;10*ones(2*nx,1)];
umina = [umin;zeros(2*nx,1)];

nxa = nx + nu; %% x is enhaced to include u-1
nua = nu + 2*nx; %% u is enhaced to include soft constraints

%% Intermediate and End Stage Hessian
[~,P] = dlqr(Aa,Ba,Q,R,M);
H = [Q, M; M', R];

%% MPC setup
N = 10;

%% FORCES multistage form - zi = [xi, ui] for i=1...N-1 and zN = xN

stages = MultistageProblem(N);

for i = 1:N
    % initial stage
    if( i==1 )
        
        % dimension
        stages(i).dims.n = nxa+nua; % number of stage variables
        stages(i).dims.r = 2*nxa;  % number of equality constraints        
        stages(i).dims.l = nua; % number of lower bounds
        stages(i).dims.u = nua; % number of upper bounds
        stages(i).dims.p = 2*nu + 2*nx;  % number of polytopic constraints
        stages(i).dims.q = 0;     % number of quadratic constraints
        
        % cost
        stages(i).cost.H = H; % blkdiag(Q,R);
        stages(i).cost.f = [zeros(nxa,1);zeros(nu,1);ones(2*nx,1)];
        
        % lower bounds
        stages(i).ineq.b.lbidx = nxa+1:stages(i).dims.n; % lower bound acts on these indices
        stages(i).ineq.b.lb = umina; % lower bound for this stage variable
        
        % upper bounds
        stages(i).ineq.b.ubidx = nxa+1:stages(i).dims.n; % upper bound acts on these indices
        stages(i).ineq.b.ub = umaxa; % upper bound for this stage variable
        
        % equality constraints
        stages(i).eq.C = [eye(nxa), zeros(nxa,nua); Aa, Ba];
        params(1) = newParam('z1',1,'eq.c'); % RHS of first eq. constr. is a parameter: [x0, 0]
        
        % inequality constraints
        stages(i).ineq.p.A = [zeros(nu,nx),eye(nu),-eye(nu),zeros(nu,nx),zeros(nu,nx);...
                              zeros(nu,nx),-eye(nu),eye(nu),zeros(nu,nx),zeros(nu,nx);...
                              C,zeros(nx,nu),zeros(nx,nu),-eye(nx),zeros(nx);...
                              -C,zeros(nx,nu),zeros(nx,nu),zeros(nx),-eye(nx)];
        stages(i).ineq.p.b = [dumax; dumax; xmax; -xmin];
        
    end
    
    % stages along horizon
    if( i>1 && i<N )       
        
        % dimension
        stages(i).dims.n = nxa+nua; % number of stage variables
        stages(i).dims.r = nxa;    % number of equality constraints        
        stages(i).dims.l = nua; % number of lower bounds
        stages(i).dims.u = nua; % number of upper bounds
        stages(i).dims.p = 2*nu + 2*nx;     % number of polytopic constraints
        stages(i).dims.q = 0;     % number of quadratic constraints
        
        % cost
        stages(i).cost.H = H; % blkdiag(Q,R);
        stages(i).cost.f = [zeros(nxa,1);zeros(nu,1);ones(2*nx,1)];
        
        % lower bounds
        stages(i).ineq.b.lbidx = 1+nxa:stages(i).dims.n; % lower bound acts on these indices
        stages(i).ineq.b.lb = umina; % lower bound for this stage variable
        
        % upper bounds
        stages(i).ineq.b.ubidx = 1+nxa:stages(i).dims.n; % upper bound acts on these indices
        stages(i).ineq.b.ub = umaxa; % upper bound for this stage variable
        
        % equality constraints
        stages(i).eq.C = [Aa, Ba];
        stages(i).eq.c = zeros(nxa,1);
        if( i==2 )
            stages(i).eq.D = [zeros(nxa,nxa+nua); -eye(nxa), zeros(nxa,nua)];
        else
            stages(i).eq.D = [-eye(nxa), zeros(nxa,nua)];
        end
        
        % inequality constraints
        stages(i).ineq.p.A = [zeros(nu,nx),eye(nu),-eye(nu),zeros(nu,nx),zeros(nu,nx);...
                              zeros(nu,nx),-eye(nu),eye(nu),zeros(nu,nx),zeros(nu,nx);...
                              C,zeros(nx,nu),zeros(nx,nu),-eye(nx),zeros(nx);...
                              -C,zeros(nx,nu),zeros(nx,nu),zeros(nx),-eye(nx)];
        stages(i).ineq.p.b = [dumax; dumax; xmax; -xmin];
        
    end
    
    % final stage
    if( i==N )
        
        % dimension
        stages(i).dims.n = nxa+2*nx;    % number of stage variables
        stages(i).dims.r = 0;     % number of equality constraints        
        stages(i).dims.l = 2*nx;    % number of lower bounds
        stages(i).dims.u = 2*nx;    % number of upper bounds
        stages(i).dims.p = 2*nx;     % number of polytopic constraints
        stages(i).dims.q = 0;     % number of quadratic constraints
        
        % cost
        stages(i).cost.H = [P,zeros(nxa,nx),zeros(nxa,nx);zeros(nx,nxa),zeros(nx),zeros(nx);zeros(nx,nxa),zeros(nx),zeros(nx)];
        stages(i).cost.f = [zeros(nxa,1);ones(2*nx,1)];
        
        % lower bounds
        stages(i).ineq.b.lbidx = nxa+1:stages(i).dims.n; % lower bound acts on these indices
        stages(i).ineq.b.lb = zeros(2*nx,1); % lower bound for this stage variable
        
        % upper bounds
        stages(i).ineq.b.ubidx = nxa+1:stages(i).dims.n; % upper bound acts on these indices
        stages(i).ineq.b.ub = 10*ones(2*nx,1); % upper bound for this stage variable
        
        % equality constraints        
        stages(i).eq.D = [-eye(nxa),zeros(nxa,nx),zeros(nxa,nx)];
        
        % inequality constraints
        stages(i).ineq.p.A = [C,zeros(nx,nu),-eye(nx),zeros(nx);...
                              -C,zeros(nx,nu),zeros(nx),-eye(nx)];
        stages(i).ineq.p.b = [xmax; -xmin];
        
    end
end

%% define outputs of the solver
outputs = newOutput('u1',1,nxa+1:nxa+nua);

%% solver settings
codeoptions = getOptions('APCV290114');

%% generate code
generateCode(stages,params,codeoptions,outputs);

%% simulate
zinit = [9.9; 0.9; 0.9; -0.1; -1; 1];
kmax = 30;
X = zeros(nx+nu,kmax+1); X(:,1) = zinit;
U = zeros(nu+2*nx,kmax);
problem.z1 = zeros(2*(nx+nu),1);

for k = 1:kmax
    problem.z1(1:nx+nu) = X(:,k);
    [solverout,exitflag,info] = APCV290114(problem);
    if( exitflag == 1 )
        U(:,k) = solverout.u1;
    else
        info
        error('Some problem in solver');
    end
    X(1:nx,k+1) = A*X(1:nx,k) + B*U(1:nu,k);
    X(nx+1:nx+nu,k+1) = U(1:nu,k);
end

%% plot
figure(1); clf;
subplot(2,1,1); grid on; title('states'); hold on;
plot([1 kmax], [xmax xmax]', 'r--'); plot([1 kmax], [xmin xmin]', 'r--');
ylim(1.1*[-5,max(xmax)]); stairs(1:kmax,X(:,1:kmax)');

subplot(2,1,2);  grid on; title('input'); hold on;
plot([1 kmax], [umax umax]', 'r--'); plot([1 kmax], [umin umin]', 'r--');
ylim(1.1*[min(umin),max(umax)]); stairs(1:kmax,U(1:nu,1:kmax)');
