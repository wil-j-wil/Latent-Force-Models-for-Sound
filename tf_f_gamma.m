% Slightly edited version of TF_F from the LFM toolbox. Now includes
% gamma parameter.

% TF_F Transcription factor model (1st order non-linear LFM) dynamics
% 
% Syntax:
%   f = tf_f(x,t,param)
%
% In:
%   x     - State vector
%   t     - Time instance
%   param - Model parameters (see below for details)    
% Out:
%   f     - Jacobian of f(x(t),t) wrt x
%
% Description:
%   Dynamics of transcriptional factor, described by the equations
%   
%   dx_i(t)/dt = B_i + \sum_{r=1}^R S_{i,r} g_i(u_r(t)) - D_i x(t),
% 
%   for outputs x_i(t), i=1,...,N and latent forces u_r(t), r=1,...,R,
%   which modelled with LTI SDE models specified by the user in param
%   structure. The non-linear function g_i(.) is also specified in param.
%

% Copyright (C) 2012 Jouni Hartikainen
%
% This software is distributed under the GNU General Public 
% Licence (version 3 or later); please refer to the file 
% Licence.txt, included with the software, for details.

% Updated in 2017 by William J. Wilkinson

function f = tf_f_gamma(x,t,param)
    
    % Extract needed parameters
    A = param.A;
    B = param.B;
    D = param.D;
    S = param.S;    
    d = param.d;
    F = param.F;
    H = param.H;
    g = param.g_func;
    g_param = param.g_param;
    
    gamma = param.gamma; % WJW - Mar 2017

    % Storage for time derivative
    f = zeros(size(x));   

    % Evaluate transformed forces
    gf = feval(g,H*x(d+1:end,:),g_param);

    % TF dynamics
    for i1 = 1:size(x,2)
        %f(1:d,i1) = B + S.*gf(:,i1) - D.*x(1:d,i1); % old version        
        f(1:d,i1) = S.*gf(:,i1) - D.*real(x(1:d,i1).^(gamma)); % WJW - Mar 2017
    end
    
    
    % GP dynamics
    f(d+1:end,:) = F*x(d+1:end,:);