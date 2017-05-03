% A slightly edited version of E_SQRT_GF from the LFM toolbox. Now includes a counter.

% E_SQRT_GF Negative log likelihood with sqrt Gaussian filter
% 
% Syntax:
%  E = E_sqrt_gf(theta,param)
%
% In:
%   theta - Parameters as dx1 vector
%   param - Parameter structure (see below for details)
%      
% Out:
%   E - Negative log likelihood
% 
% Description:
%  Calculates the negative log likelihood of data with a partially observed
%  SDE model 
% 
%  dx(t)/dt = f(x(t),t)   + L(x(t),t) w(t)
%       y_k = h_k(x(t_k)) + r_k, r_k ~ N(0,R_k), k = 1,...,T
% 
%  using a square-root form continuous-discrete Gaussian filter. 
%
% Copyright (C) 2012 Jouni Hartikainen, Simo Särkkä
%
% This software is distributed under the GNU General Public 
% Licence (version 3 or later); please refer to the file 
% Licence.txt, included with the software, for details.

% Updated in 2017 by William J. Wilkinson

function [E] = E_sqrt_gf_count(theta, param)

Y          = param.Y;           % Measurements
mind       = param.mind;        % Measurement index vector
f_func     = param.f_func;      % Dynamic model function
f_param    = param.f_param;     % Parameters of f 
q_func     = param.q_func;      % Dispersion function
q_param    = param.q_param;     % Parameters of q   
param_func = param.param_func;  % Function for setting parameter vectors w. theta
int_method = param.int_method;  % Moment integration method
h_func     = param.h_func;      % Measurement model function
h_param    = param.h_param;     % Parameters of h
R          = param.R;           % Measurement noise covariance (matrix or function returning matrix)
dt         = param.dt;          % Time step between measurements   
t_min      = param.t_min;       % Time at initial time step
ltr_ind    = param.ltr_ind;     % Indicator vector for log-transformed parameters

% Parameters of R if it's a function
if isfield(param,'n_param')
    n_param = param.n_param;
else
    n_param = [];
end

if ~isnumeric(R)
    R = R(theta,n_param);
end

% Function used to solve ODEs of m and P wrt time
if isfield(param,'ode_func')
    ode_func = param.ode_func;
else
    % This is usually the best choice
    ode_func = @ode45;
    %ode_func = @ode15s; % WJW - Mar 2017
end

% Parameters of ODE solver
if isfield(param,'ode_opt')
    opt = param.ode_opt;
else
    opt = odeset;
end
%opt.MaxOrder = 5; % Could set to 1 to make faster. Not sure how it would affect accuracy. WJW - Mar 2017

if isfield(param,'start_ind')
    start_ind = param.start_ind;
else
    start_ind = 0;
end

% Parameters of moment integration function
if isfield(param,'int_param')
    int_param = param.int_param;
else
    int_param = {};
end

% Exponentiate the log transformed parameters
theta(ltr_ind) = exp(theta(ltr_ind));

% Add prior contribution
E_prior = 0;
if isfield(param,'e_prior')
    for i = 1:length(theta)
        E_prior = E_prior - real(feval(param.e_prior{i}.fh.lp,theta(i),param.e_prior{i}));   
    end
end

% Set theta to param, f_param and q_param
[param, f_param, q_param] = feval(param_func,param,f_param,q_param,theta);

% Dispersion matrix
Q = feval(q_func,[],[],q_param);

% This notation is used with Gaussian filters
Q = Q*Q';

% Prior mean and covariance
M0 = param.M0;
P0 = param.P0;

m = M0;
P = P0;
A = chol(P)';

n = size(m,1);
steps = length(mind);

% Space for conditional measurement likelihoods
EE = zeros(1,steps);

% Form X and W for the chosen sigma point rule
[XI,WM,WC] = feval(int_method,n,int_param);

% Parameters for the moment derivative function to be integrated
ode_param = {XI,WM,WC,f_func,f_param,Q,t_min};

% Measurement counter
mc = 1;
for k=1:steps
    %k
    mA = [m A];
    if k > start_ind
        vmA = mA2vec(mA,n);
        t = t_min+dt*(k-1);
        Tspan = [t t+dt];
        [sol,y] = feval(ode_func,@(t,y) myode(t,y,@dmA_sqrt,ode_param),Tspan,vmA,opt);
        vmA = y(end,:)';
        [mA,m,A] = vec2mA(vmA,n);
    end
    
    if mind(k) == 1
        % Linear update
        if isnumeric(h_func)
            [m,P,A,W,LLH] = srkf_update(m,A,Y(:,mc),h_func,chol(R)');
        % Non-linear update
        else
            [m,P,A,W,LLH] = sckf_update(m,A,Y(:,mc),h_func,R,h_param);
        end
        mc = mc + 1;
        EE(k) = LLH;
    end
    %k
end
E = -sum(EE) + E_prior; 
global count
count = count + 1