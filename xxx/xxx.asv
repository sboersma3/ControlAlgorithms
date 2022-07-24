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
[r,y,x,u,d] = DefineSignals(ops);

% time-loop
CPUTime = zeros(1,ops.N);
for kk=1:ops.N
    
    tic

    % propagate system one time step ahead
    x(:,kk+1)       = f(x(:,kk), u(:,kk), d(:,kk), p, ops.h);
    y(:,kk)         = g(x(:,kk), u(:,kk), d(:,kk), p, ops.h);
    
    % controller
    [u(:,kk+1),ops] = Controller(r(:,kk),x(:,kk),u(:,kk),d(:,kk),p,ops);

    CPUTime(kk)     = toc;

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


ops.nr = 1;       % #references
ops.ny = 1;       % #measurements
ops.nx = 1;       % #states
ops.nu = 1;       % #control signals
ops.nd = 1;       % #disturbances

end

function ops = DefineControllerSettings(p,ops)

% controller parameters


end

function [r,y,x,u,d] = DefineSignals(ops)

r          = zeros(ops.nr,ops.N);                                       % reference
y          = zeros(ops.ny,ops.N);                                       % measurement
x          = zeros(ops.nx,ops.N+1);                                     % state
x(:,1)     = zeros(ops.nx,1);                                           % initial state
u          = zeros(ops.nu,ops.N+1);                                     % control signal
u(:,1)     = zeros(ops.nu,1);                                           % initial control signal
d          = zeros(ops.nd,ops.N);                                       % disturbance

end

function [u,ops] = Controller(r,x,u,d,p,ops)


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

ki = [];

end

function y = ghat(x,u,d,p,h)

y = [];

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

ki = [].*(1+0.05*randn(size(x,1),1));

end

function y = g(x,u,d,p,h)

y = [];
y = y.*(1+0.01*randn(size(y,1),1));

end