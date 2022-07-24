clear;clc

% NMPC example using Casadi and the single shooting method

addpath(genpath('C:\casadi-windows-matlabR2016a-v3.5.5'));
addpath('bin');

% model parameters
[p,ops]     = DefineModelParameters();

% simulation parameters
ops.N       = 150;                      % #samples in simulation
ops.h       = .1;                       % sample period
ops.t       = 0:ops.h:ops.N*ops.h;      % time

% controller settings
ops         = DefineControllerSettings(p,ops);

% signals
[r,y,x,u,d,xhorizon,uopt] = DefineSignals(ops);

% time-loop
CPUTime = zeros(1,ops.N);
for kk=1:ops.N
    
    tic
    % propagate system one time step ahead
    x(:,kk+1)  = f(x(:,kk), u(:,kk), d(:,kk), p, ops.h);
    y(:,kk)    = g(x(:,kk), u(:,kk), d(:,kk), p, ops.h);
    
    [u(:,kk+1),uopt(:,:,kk+1),xhorizon(:,:,kk),ops] = ...
        Controller(r(:,kk:kk+ops.k1-1),x(:,kk),u(:,kk),d(:,kk:kk+ops.k1-1),p,ops);
    
    CPUTime(kk)  = toc;
end

disp(' ')
disp(['mean CPU time: ',num2str(mean(CPUTime)),' (sec).'])
disp(' ')

PlotResults

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% help functions
function [p,ops] = DefineModelParameters()

p.wn   = 1.23;    % model parameter 1
p.zeta = 5;       % model parameter 2
ops.nr = 1;       % #references
ops.ny = 1;       % #measurements
ops.nx = 2;       % #states
ops.nu = 1;       % #control signals
ops.nd = 0;       % #disturbances

end

function ops = DefineControllerSettings(p,ops)

% controller parameters
ops.k1     = 50;                    % #samples in the prediction horizon
ops.u_min  = -3/2;                  % lower bound on controllable inputs
ops.u_max  = 3/2;                   % upper bound on controllable inputs
ops.x_min  = [-6/2;-6/2];           % lower bound on states    
ops.x_max  = [6/2;6/2];             % upper bound on states 
ops.du_max = .10;                   % lower and upper bound change of u 
ops.u0     = zeros(ops.nu,ops.k1);  % initial conditions decision variables

% MPC solver options
ops.opts.ipopt.print_level              = 3;
ops.opts.ipopt.max_iter                 = 1000;
ops.opts.ipopt.nlp_scaling_method       = 'gradient-based';
ops.opts.ipopt.warm_start_init_point    = 'yes';


% cost and nonlinear constraints
ops        = DefineCostConstraints(p,ops);


end

function ops = DefineCostConstraints(p,ops)

% us = [ u1_0 u1_1 .. u1_k1-1 ; 
%        u2_0 u2_1 .. u2_k1-1 ; 
%        u3_0 u3_1 .. u3_k1-1 ...], with ui_j = ui(j)

% Us = [u1(0) u1(1) .. u1(k1-1) u2(0) u2(1) ... u2(k1-1)]

x0  = casadi.SX.sym('x0',ops.nx,1);         % initial condition
us  = casadi.SX.sym('us',ops.nu,ops.k1);    % control signal
dus = diff(us,1,2);                         % difference control signal
Us  = us'; Us = Us(:);                      % vectorized control signal Us
ds  = casadi.SX.sym('ds',ops.nd,ops.k1);    % disturbance signal
xs  = casadi.SX.sym('xs',ops.nx,ops.k1+1);  % state
ys  = casadi.SX.sym('ys',ops.ny,ops.k1);    % measurement
rs  = casadi.SX.sym('rs',ops.nr,ops.k1);    % reference

% propagate the model symbolically forward in time over the prediction horizon
xs(:,1) = x0;
for ll=1:ops.k1
    xs(:,ll+1) = fhat(xs(:,ll), us(:,ll), ds(:,ll), p, ops.h);
    ys(:,ll)   = ghat(xs(:,ll), us(:,ll), ds(:,ll), p, ops.h);
end

% choose here a cost Js, 
% (non)linear inequality contraints cs(x) < 0 and 
% (non)linear equality contraints ceqs(x) = 0
Js = 0;
for ll=1:ops.k1
    Js          = Js + (rs(:,ll)-ys(:,ll))'*20*eye(ops.ny)*(rs(:,ll)-ys(:,ll)) + ...
                                    us(:,ll)'*.1*eye(ops.ny)*us(:,ll);
end
cs              = [xs(1,2:end)-ops.x_max(1)...
                    xs(2,2:end)-ops.x_max(2)...
                    -xs(1,2:end)+ops.x_min(1)...
                    -xs(2,2:end)+ops.x_min(2)];
ceqs            = [];                     

ops.xs          = casadi.Function('xs',{x0,us,ds},{xs});
ops.F           = casadi.Function('F',{x0,Us,ds,rs},{Js,cs,ceqs}, {'x0','Us','ds','rs'}, {'J', 'c','ceq'});

ops.us          = us;
ops.Us          = Us;
ops.dus         = dus;

% bounds on the control signals
ops.lbu = [];
ops.ubu = [];
for ll=1:ops.nu
    ops.lbu = [ops.lbu;ops.u_min(ll)*ones(ops.k1,1)];
    ops.ubu = [ops.ubu;ops.u_max(ll)*ones(ops.k1,1)];
end

ops.lbg = [-inf*ones(size(cs,2),1) ; zeros(size(ceqs,2),1)];
ops.ubg = [zeros(size(cs,2),1) ; zeros(size(ceqs,2),1)];

for ll=1:ops.nu
    ops.lbg   = [ops.lbg ; -ops.du_max(ll)*ones(ops.k1,1)];
    ops.ubg   = [ops.ubg ; ops.du_max(ll)*ones(ops.k1,1)];
end

end

function [r,y,x,u,d,xhorizon,uopt] = DefineSignals(ops)

r          = [1.0*ones(1,floor((ops.N+ops.k1)/2)) ...
                1.3*ones(1,ceil((ops.N+ops.k1)/2))];                    % reference
y          = zeros(ops.ny,ops.N);                                       % measurement
x          = zeros(ops.nx,ops.N+1);                                     % state
x(:,1)     = zeros(ops.nx,1);                                           % initial state
u          = zeros(ops.nu,ops.N+1);                                     % control signal
u(:,1)     = ops.u0(:,1);                                               % initial control signal
d          = zeros(ops.nd,ops.N+ops.k1);                                % disturbance

xhorizon   = zeros(ops.nx,ops.k1+1,ops.N);                              % state predictions in the horizon
uopt       = zeros(ops.nu,ops.k1,ops.N);                                % control signals in the horizon

end

function [u,uopt,xhorizon,ops] = Controller(r,x,u,d,p,ops)

% update linear inequality constraint
ops.dUs           = horzcat(ops.us(:,1)-u,ops.dus)'; ops.dUs = ops.dUs(:);

% solve the optimization problem
[uopt,res.output] = Nmpc(r,x,d,p,ops);

u                 = uopt(:,1);

ops.lab_g0        = res.output.lam_g;        % for warmstart
ops.lab_x0        = res.output.lam_x;        % for warmstart
ops.u0            = uopt;                    % for warmstart

% state evolutions in the prediction horizon
xhorizon          = full(ops.xs(x,uopt,d));

end

%% controller model
function dx = fhat(x,u,d,p,h)

k1  = Fhat(x,u,d,p);
k2  = Fhat(x + h/2 * k1,u,d,p);
k3  = Fhat(x + h/2 * k2,u,d,p);
k4  = Fhat(x + h * k3,u,d,p);
dx  = x + h/6*(k1 + 2*k2 + 2*k3 + k4);

end

function ki = Fhat(x,u,d,p)

ki = [x(2);-2*p.zeta*p.wn*x(2)^2-p.wn^2*x(1)] + [0;p.wn^2*u(1)];

end

function y = ghat(x,u,d,p,h)

y = x(1);

end


%% system
function dx = f(x,u,d,p,h)

k1  = F(x,u,d,p);
k2  = F(x + h/2 * k1,u,d,p);
k3  = F(x + h/2 * k2,u,d,p);
k4  = F(x + h * k3,u,d,p);
dx  = x + h/6*(k1 + 2*k2 + 2*k3 + k4);

end

function ki = F(x,u,d,p)

ki = [x(2);-2*p.zeta*p.wn*x(2)^2-p.wn^2*x(1)].*(1+0.05*randn(size(x,1),1)) + [0;p.wn^2*u(1)];

end

function y = g(x,u,d,p,h)

y = x(1);
y = y.*(1+0.01*randn(size(y,1),1));

end