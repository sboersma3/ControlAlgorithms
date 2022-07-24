clear;clc

% NMPC example using Casadi and the single shooting method
% The external disturbance (d) is not modelled in the NMPC

addpath(genpath('C:\casadi-windows-matlabR2016a-v3.5.5'));
addpath('bin');

% MPC solver options
ops.opts.ipopt.print_level              = 3;
ops.opts.ipopt.max_iter                 = 1000;
ops.opts.ipopt.nlp_scaling_method       = 'gradient-based';
ops.opts.ipopt.warm_start_init_point    = 'yes';

% model parameters
p.wn   = 1.23;    % model parameter 1
p.zeta = 5;       % model parameter 2
ops.nx = 2;       % #states
ops.ny = 1;       % #measurements
ops.nu = 1;       % #control signals
ops.nr = 1;       % #references
ops.nd = 1;       % #disturbances

% simulation parameters
ops.N   = 150;                      % #samples in simulation
ops.h   = .1;                       % sample period
ops.t   = 0:ops.h:ops.N*ops.h;      % time

% controller parameters
ops.Np     = 20;                    % #samples in the prediction horizon
ops.u_min  = -3/2;                  % lower bound on controllable inputs
ops.u_max  = 3/2;                   % upper bound on controllable inputs
ops.x_min  = [-6/2;-6/2];           % lower bound on states    
ops.x_max  = [6/2;6/2];             % upper bound on states 
ops.du_max = [.10 .10];             % lower and upper bound change of u 
ops.u0     = zeros(ops.nu,ops.Np);  % initial conditions decision variables

% cost and nonlinear constraints
ops        = CostConstraintFunctions(p,ops);

% signals
r          = [1.0*ones(1,floor((ops.N+ops.Np)/2)) ...
                1.0*ones(1,ceil((ops.N+ops.Np)/2))];                    % reference
d          = [zeros(1,floor((ops.N+ops.Np)/3)) ...
                -.7*ones(1,ceil(2*(ops.N+ops.Np)/3))];                  % disturbance
d          = lsim(tf(1,[.5 1]),d,0:ops.h:(ops.N+ops.Np)*ops.h-ops.h)';  % filter disturbance
x          = zeros(ops.nx,ops.N+1);                                     % state
y          = zeros(ops.ny,ops.N);                                       % measurement
u          = zeros(ops.nu,ops.N+1);                                     % control signal
u(:,1)     = ops.u0(:,1);                                               % initial control signal
xhorizon   = zeros(ops.nx,ops.Np+1,ops.N);                              % state predictions in the horizon
yhorizon   = zeros(ops.ny,ops.Np,ops.N);                                % measurement prediction in the horizon
uopt       = zeros(ops.nu,ops.Np,ops.N);                                % control signals in the horizon

res.output = cell(1,ops.N);
CPUTime    = zeros(1,ops.N);

%%%%% time-loop %%%%%
for kk=1:ops.N
    
    tic
    % propagate model one time step ahead
    x(:,kk+1)  = f(x(:,kk), u(:,kk), d(:,kk), p, ops.h);
    y(:,kk)    = g(x(:,kk), u(:,kk), d(:,kk), p, ops.h);
    
    % update linear inequality constraint
    ops.dUs = horzcat(ops.us(:,1)-u(:,kk),ops.dus)'; ops.dUs = ops.dUs(:);
    
    % solve the optimization problem
    [uopt(:,:,kk+1),res.output{kk}] = ...
        Nmpc(x(:,kk),0*d(:,kk:kk+ops.Np-1),r(:,kk:kk+ops.Np-1),p,ops);
        
    u(:,kk+1)    = uopt(:,1,kk+1);
    
    ops.lab_g0   = res.output{kk}.lam_g;        % for warmstart
    ops.lab_x0   = res.output{kk}.lam_x;        % for warmstart
    ops.u0       = uopt(:,:,kk+1);              % for warmstart
    
    % state evolutions in the prediction horizon
    xhorizon(:,:,kk) = full(ops.xs(x(:,kk),uopt(:,:,kk+1),d(:,kk:kk+ops.Np-1)));
    yhorizon(:,:,kk) = full(ops.ys(x(:,kk),uopt(:,:,kk+1),d(:,kk:kk+ops.Np-1)));
    
    CPUTime(kk)  = toc;
end

disp(' ')
disp(['mean CPU time: ',num2str(mean(CPUTime)),' (sec).'])
disp(' ')

% plot
figure(1);clf;
subplot(3,2,1)
stairs(ops.t(1,1:kk),x(1,1:kk),'linewidth',1.5);grid;hold on
stairs(ops.t(1,1:kk),ops.x_max(1)*ones(1,kk),'r--');
stairs(ops.t(1,1:kk),ops.x_min(1)*ones(1,kk),'r--');
ylabel('$x_1$','fontsize',18,'interpreter','latex')
subplot(3,2,2)
stairs(ops.t(1,1:kk),x(2,1:kk),'linewidth',1.5);grid;hold on
stairs(ops.t(1,1:kk),ops.x_max(2)*ones(1,kk),'r--');
stairs(ops.t(1,1:kk),ops.x_min(2)*ones(1,kk),'r--');
ylabel('$x_2$','fontsize',18,'interpreter','latex')
subplot(3,2,3)
stairs(ops.t(1,1:kk),u(1,1:kk),'linewidth',1.5);grid;hold on
stairs(ops.t(1,1:kk),ops.u_max(1)*ones(1,kk),'r--');
stairs(ops.t(1,1:kk),ops.u_min(1)*ones(1,kk),'r--');
ylabel('$u_1$','fontsize',18,'interpreter','latex')
xlabel('$t$ (s)','fontsize',18,'interpreter','latex')
subplot(3,2,4)
stairs(ops.t(1,1:kk),d(1,1:kk),'linewidth',1.5);grid
ylabel('$d_1$','fontsize',18,'interpreter','latex')
subplot(3,2,6)
stairs(ops.t(1,1:kk),y(1,1:kk),'linewidth',1.5);grid;hold on
stairs(ops.t(1,1:kk),r(1,1:kk),'k--','linewidth',1.5);
ylabel('$y_1$','fontsize',18,'interpreter','latex')
xlabel('$t$ (s)','fontsize',18,'interpreter','latex')
%
function ops = CostConstraintFunctions(p,ops)

% us = [ u1_0 u1_1 .. u1_Np-1 ; 
%        u2_0 u2_1 .. u2_Np-1 ; 
%        u3_0 u3_1 .. u3_Np-1 ...], with ui_j = ui(j)

% Us = [u1(0) u1(1) .. u1(Np-1) u2(0) u2(1) ... u2(Np-1)]

x0  = casadi.SX.sym('x0',ops.nx,1);         % initial condition
us  = casadi.SX.sym('us',ops.nu,ops.Np);    % control signal
dus = diff(us,1,2);                         % difference control signal
Us  = us'; Us = Us(:);                      % vectorized control signal Us
ds  = casadi.SX.sym('ds',ops.nd,ops.Np);    % disturbance signal
xs  = casadi.SX.sym('xs',ops.nx,ops.Np+1);  % state
ys  = casadi.SX.sym('ys',ops.ny,ops.Np);    % measurement
rs  = casadi.SX.sym('rs',ops.nr,ops.Np);    % reference

% propagate the model symbolically forward in time over the prediction horizon
xs(:,1) = x0;
for ll=1:ops.Np
    xs(:,ll+1) = fhat(xs(:,ll), us(:,ll), ds(:,ll), p, ops.h);
    ys(:,ll)   = ghat(xs(:,ll), us(:,ll), ds(:,ll), p, ops.h);
end

% choose here a cost Js, 
% (non)linear inequality contraints cs < 0 and 
% (non)linear equality contraints ceqs = 0
Js = 0;
for ll=1:ops.Np
    Js          = Js + (rs(:,ll)-ys(:,ll))'*20*eye(ops.ny)*(rs(:,ll)-ys(:,ll)) + ...
                                    us(:,ll)'*.1*eye(ops.ny)*us(:,ll);
end
cs              = [xs(1,2:end)-ops.x_max(1)...
                    xs(2,2:end)-ops.x_max(2)...
                    -xs(1,2:end)+ops.x_min(1)...
                    -xs(2,2:end)+ops.x_min(2)];
ceqs            = [];                     

ops.xs          = casadi.Function('xs',{x0,us,ds},{xs});
ops.ys          = casadi.Function('ys',{x0,us,ds},{ys});
ops.F           = casadi.Function('F',{x0,Us,ds,rs},{Js,cs,ceqs}, {'x0','Us','ds','rs'}, {'J', 'c','ceq'});

ops.us          = us;
ops.Us          = Us;
ops.dus         = dus;

% bounds on the control signals
ops.lbu = [];
ops.ubu = [];
for ll=1:ops.nu
    ops.lbu = [ops.lbu;ops.u_min(ll)*ones(ops.Np,1)];
    ops.ubu = [ops.ubu;ops.u_max(ll)*ones(ops.Np,1)];
end

ops.lbg = [-inf*ones(size(cs,2),1) ; zeros(size(ceqs,2),1)];
ops.ubg = [zeros(size(cs,2),1) ; zeros(size(ceqs,2),1)];

for ll=1:ops.nu
    ops.lbg   = [ops.lbg ; -ops.du_max(ll)*ones(ops.Np,1)];
    ops.ubg   = [ops.ubg ; ops.du_max(ll)*ones(ops.Np,1)];
end

end

%% contorller model
function dx = fhat(x,u,d,p,h)

k1  = Fhat(x,u,d,p,h);
k2  = Fhat(x + h/2 * k1,u,d,p,h);
k3  = Fhat(x + h/2 * k2,u,d,p,h);
k4  = Fhat(x + h * k3,u,d,p,h);
dx  = x + h/6*(k1 + 2*k2 + 2*k3 + k4);

end

function ki = Fhat(x,u,d,p,h)

ki = [x(2);-2*p.zeta*p.wn*x(2)^2-p.wn^2*x(1)] + [0;p.wn^2*u(1)] + [0;d(1)];

end

function y = ghat(x,u,d,p,h)

y = x(1);

end


%% system
function dx = f(x,u,d,p,h)

k1  = F(x,u,d,p,h);
k2  = F(x + h/2 * k1,u,d,p,h);
k3  = F(x + h/2 * k2,u,d,p,h);
k4  = F(x + h * k3,u,d,p,h);
dx  = x + h/6*(k1 + 2*k2 + 2*k3 + k4);

end

function ki = F(x,u,d,p,h)

ki = [x(2);-2*p.zeta*p.wn*x(2)^2-p.wn^2*x(1)].*(1+0.05*randn(size(x,1),1)) + [0;p.wn^2*u(1)] + [0;d(1)];

end

function y = g(x,u,d,p,h)

y = x(1)*(1+0.01*randn);

end