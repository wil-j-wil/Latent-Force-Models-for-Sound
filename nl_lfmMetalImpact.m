%% Analysis Data
[Actual,fs] = audioread('metal-impact.wav'); % the audio data
[allMag, allFreq] = readSpearFile('metal-impact.txt'); % sinusoidal analysis data from Spear
hopSize     = fs / 100; % Spear exports data with 10ms frames
saveResults = 1;    % save the results to a mat file?
matNameMag  = 'resultsMetalImpact_2.mat'; % name of the .mat file saved if saveResults = 1

%%%%%%%%%%%%%
% when running the model, if the count stops increasing then optimisation
% is likely stuck. Try a different choice of gamma, fewer modes, different
% intial parameter guess, or change the noise assumptions (line 235/237)
%%%%%%%%%%%%%
 
% first we must select our modes of vibration, and format the data
 
 % pick the _ modes of vibration
 numModes = 8;
 
 % peak picking polynomial order (change this if the wrong modes are being selected)
 polyOrder = 12;
 
 % sample every _ data points (to reduce data size)
 sfactor = 3;
 
 % "linearity" measure [1/2,1] (1/2=linear decay, 1=exponential decay)
 gamma = 3/4;
 
 % maximum number of function evaluations during optimisation (the "count" variable)
 MaxFunEvals = 1500; % even when count > MaxFunEvals, the optimiser will wait until the end of the current iteration
 
 param_guess = 1; % make an inital parameter guess? 1 = use D_ and S_, 0 = try random
    % these guesses are fairly good, so the optimiser probably wont change them much
    D_ = [0.4; 1.3; 0.8; 0.8; 0.8; 0.4; 0.4; 1.3] + 0.05*rand(numModes,1); % guess for damping (partly randomised)
    S_ = [0.6; 0.6; 0.3; 0.3; 0.3; 0.3; 0.3; 0.3] + 0.05*rand(numModes,1); % guess for stiffness (partly randomised)
    
 priorityInference = 1; % 1 = smaller noise assumption for modes with larger amplitude, so they get priority
 
 
 %% Find the most prominent modes and format / scale the data
 
 % normalise so that maximum value is 0.3
 allMagMax = max(allMag');
 maxMax = max(allMagMax);
 allMag = 0.3 * allMag / maxMax;

 % plot the sinusoidal data
 figure(1); clf
 subplot(211)
 plot(allMag')
 title('Magnitude Analysis')
 xlabel('Frame Number')
 ylabel('Magnitude')
 subplot(212)
 plot(allFreq')
 title('Frequency Analysis')
 xlabel('Frame Number')
 ylabel('Frequency (Hz)')
 
 % calculate the median frequency values, since we will be using fixed frequencies
 allFreq(find(allFreq==0)) = NaN; % set zeros to NaN so that they aren't included in median calc
 medFreq = zeros(size(allMag,1),1);
 for i=1:size(allMag,1)
     medFreq(i) = nanmedian(allFreq(i,:)); % median frequency value
 end
 notNaN = find(isnan(medFreq)==0);
 allFreq = allFreq(notNaN,:); % discard empty partials
 allMag = allMag(notNaN,:);
 medFreq = medFreq(notNaN,:);

 [Bs,I] = sort(medFreq,'ascend'); % sort by frequency
 allMag = allMag(I,:);
 allFreq = allFreq(I,:);
 medFreq = medFreq(I,:);
 
 allMagMax = max(allMag'); % maximum of partial
 allMagSum = sum(allMag'); % sum of partial
 
 % create vector of evenly spaced frequency bins, and sum the sinusoids in each bin
 numBins  = size(allMag,1);
 maxFreq  = medFreq(end);
 binInd   = linspace(0,maxFreq,numBins);
 totalSum = 0;
 binSum   = zeros(numBins,1);
 for i=2:numBins
     sum_      = sum(allMagSum(find(medFreq<=binInd(i))));
     binSum(i) = sum_ - totalSum;
     totalSum  = sum_;
     if binSum(i) == 0
         binSum(i) = binSum(i-1); % we don't want any zeros
     end
 end
 binSum(1) = binSum(2);
 binSum(find(binSum==0))=0.001;% to catch zeros in the first bins
 
 % fit a polynomial to our frequency data
 warning('off','all');
  sinPoly2 = polyfit(binInd,mag2db(binSum)',polyOrder); % poly of order 10 or 12 works well
 warning('on','all');
 sinPolyDraw2 = polyval(sinPoly2, binInd); % this will become our reference
 
 sinPolyMag = db2mag(sinPolyDraw2); % adjustment in the magnitude domain
 
 % alter the magnitudes, based on sinPolyMag - i.e. flatten the frequency
 % domain content so that we can pick our modes
 allMagSumAlt = zeros(size(allMagSum));
 augMag = allMag;
 augFreq = allFreq;
 for j=2:numBins
     altInd = find(medFreq<=binInd(j) & medFreq>binInd(j-1));
     allMagSumAlt(altInd) = mag2db(allMagSum(altInd)) - sinPolyDraw2(j);
     augMag(altInd,:) = allMag(altInd,:) ./ sinPolyMag(j);
 end
 
 % "filter" the data
 binFilt = mag2db(binSum)-sinPolyDraw2';
 allMagFilt = mag2db(allMagSum)-sinPolyDraw2;
 
 % plot for comparison
 figure(100); clf
 subplot(221)
 plot(binInd,mag2db(binSum))
 hold on
 plot(binInd,sinPolyDraw2,'m--')
 xlabel('Frequency (Hz)')
 ylabel('Power/Frequency (dB/Hz)')
 title('Sum of sinusoids in evenly spaced freq. bins')
 subplot(223)
 plot(binInd,binFilt)
 xlabel('Frequency (Hz)')
 ylabel('Power/Frequency (dB/Hz)')
 title('"High-pass filtered" data, based on polynomial fit')
 subplot(222)
 plot(medFreq,mag2db(allMagSum))
 xlabel('Frequency (Hz)')
 ylabel('Power/Frequency (dB/Hz)')
 title('The actual sinusoidal data')
 subplot(224)
 plot(medFreq,allMagSumAlt)
 xlabel('Frequency (Hz)')
 ylabel('Power/Frequency (dB/Hz)')
 title('The "filtered" actual sinusoidal data')
 
 % calculate max and sum of our new, augmented data
 augMagMax = max(augMag');
 augMagSum = sum(augMag');
 
 % sort by sum
 [Bs,I] = sort(augMagSum,'descend');
 augMag = augMag(I,:);
 augFreq = augFreq(I,:);
 allMag = allMag(I,:);
 allFreq = allFreq(I,:);
 medFreq = medFreq(I,:);
 
 residualMag = allMag(numModes+1:end,:);
 residualFreq = allFreq(numModes+1:end,:);
 residualMedFreq = medFreq(numModes+1:end,:);
 modeMag = allMag(1:numModes,:); % just keep the _ largest partials
 modeFreq = allFreq(1:numModes,:);
 modeMedFreq = medFreq(1:numModes,:);
 
 % sort by frequency again
 [Bs,I] = sort(modeMedFreq,'ascend');
 modeMag = modeMag(I,:);
 modeFreq = modeFreq(I,:);
 modeMedFreq = modeMedFreq(I,:);
 [Bs,I] = sort(residualMedFreq,'ascend');
 residualMag = residualMag(I,:);
 residualFreq = residualFreq(I,:);
 residualMedFreq = residualMedFreq(I,:);
 
 % plot again
 figure(2); clf
 subplot(211)
 plot(modeMag')
 title('Selected Modes - Magnitude')
 xlabel('Frame Number')
 ylabel('Magnitude')
 subplot(212)
 plot(modeFreq')
 title('Selected Modes - Frequency')
 xlabel('Frame Number')
 ylabel('Frequency (Hz)')
 
 figure(3); clf
 plot(residualMag')
 title('Residual')
 %inSig = spearMag;
 
 % re-scaling so that max value of all modes is 0.3
 inSig = zeros(size(modeMag));
 maxOfMode = zeros(1,size(modeMag,1));
 scaleFactor = zeros(1,size(modeMag,1));
 for i=1:numModes
    maxOfMode(i) = max(modeMag(i,:));
    inSig(i,:) = modeMag(i,:) .* 0.3 ./ maxOfMode(i);
    scaleFactor(i) = 0.3 ./ maxOfMode(i);
 end
 
 figure(4); clf
 pIn=plot(inSig');
 title('Scaled observations')
 legnames = cell(1);
 for i=1:numModes
    legNames{i} = strcat('mode',num2str(i));
 end
 legend(pIn,legNames{1:numModes})
 
 
 
%% 1st order non-linear latent force inference with nonlinear Kalman filtering and smoothing.
% Adapted from Hartikeinen and Sarkka's LFM toolbox: http://becs.aalto.fi/en/research/bayes/lfm/
% 
% 
% The model for i = 1,...,N outputs is
% 
% dx_i(t)/dt = B_i + S_{i} g_i(u(t)) - D_i x(t)^(gamma),
% 
% where the latent force u(t) has LTI SDE prior (augmented as a part
% of the joint state-space model). 

steps = size(inSig,2);

t_min = 1;
t_max = steps;
tgrid = linspace(t_min,t_max,steps);

dt = 0.1;

% Number of outputs
d = size(inSig,1);

% Measurement noise
if priorityInference == 1
    n_sigma2 = 0.0001 ./ maxOfMode.^(4/3) % priority inference. WJW - Feb 2017
else
    n_sigma2 = 0.01^2*ones(d,1)
end

R = diag(n_sigma2);

dR = diag(R);
iR = (1./diag(R));

% Measurement indices
msteps = 1:sfactor:steps;
%msteps = 1:steps;
mind = zeros(1,steps);
mind(msteps) = 1;

% Time between measurements
dty = dt*sfactor;
mindy = ones(1,length(msteps));

% non-linearity
g_param = .1;
% "softplus" rectification function:
g_func = @(f,g_param) log(1+exp(f));
g_df   = @(f,g_param) exp(f)./(1+exp(f));
g_df2  = @(f,g_param) exp(f)./((1+exp(f)).^2);


% Inference with Gaussian filtering and smoothing:
%% First simulate the model


% Prior variance of the outputs
x_sigma2 = .0001^2;

% General parameters
param = struct;

% Parameters of f
f_func = @tf_f_gamma;
f_param = struct;
f_param.g_func  = g_func;
f_param.g_param = g_param;
f_param.g_df    = g_df;
f_param.x_sigma2 = x_sigma2;
f_param.d = d;

f_param.gamma = gamma; % WJW - Mar 2017


% LFM parameters

% ki = 1: Matern model for latent forces (in the original toolbox there is the option to pick the resonator model)
ki = 1;
f_param.ki = ki;
% Matern parameters
p = 2;
f_param.p = p;
lengthScale = 1;
sigma = 1.5;
    
theta_lf = [lengthScale;sigma];

% Initial guess for parameters
w0_lf = zeros(size(theta_lf));    

% Random ODE parameters
%Aj = .2*rand(d,1)-.1;
Aj = inSig(:,1); % initial conditions
%Bj = .1*rand(d,1);
Bj = zeros(d,1);
Dj = 0.2*rand(d,1) + 0.6; % overwritten if param_guess = 1
Sj = 0.2*rand(d,1) + 0.5; % overwritten if param_guess = 1

% use some pre-calculated params:
if param_guess == 1
    Dj = D_;
    Sj = S_;
end

theta = [Aj(:);Bj(:);Dj(:);Sj(:)];

w0    = [log(Aj);...
         log(Bj+0.000000001);...
         log(Dj);...
         log(Sj)];

w0 = [w0_lf;w0];
theta = [theta_lf;theta];


% Parameter prior
theta_prior = cell(1,length(w0));
for i = 1:length(theta)
    theta_prior{i} = prior_t;
    %theta_prior{i}.s2 = .3;
end

ltr_ind = [1:length(theta)]';

% Parameters of q
q_func = @tf_q;
q_param = struct;
q_param.d = d;

% Setup structures with theta
[param,f_param,q_param] = tf_set_params(param,f_param,q_param,theta);
M0 = param.M0;
P0 = param.P0;

x0 = gauss_rnd(M0,P0);
%x0(1:d) = Aj+Bj./Dj;
x0(1:d) = Aj;

% Dispersion matrix
[Q] = feval(q_func,[],[],q_param);

% This notation is used with Gaussian filters (with SMC above Q is used)
Q = Q*Q';

% Measurement model
h_func = [eye(d) zeros(d,size(M0,1)-d)];
h_param = [];

% Data simulation
x = x0;
X = zeros(size(x,1),steps);
T = zeros(1,steps);
Y = zeros(d,length(msteps));
t = 0;
obs_ind = 0;

mi_steps = 100;
ddt = dt/mi_steps;
XX = zeros(size(x,1),steps*mi_steps);
eu_tgrid = t_min:ddt:ddt*(steps*mi_steps-1);

for k=1:steps
    if k > 1
        for i=1:mi_steps
            t = t + ddt;
            dw = sqrt(ddt*Q) * randn(size(Q,1),1);
            dx = feval(f_func,x,t,f_param)*ddt + dw;
            x = x + dx;
            XX(:,i+(k-1)*mi_steps) = x;
        end
    end
    X(:,k) = x;
    T(:,k) = t;
    
    if ismember(k,msteps)
        obs_ind = obs_ind + 1;
        y = gauss_rnd(h_func*x,R);
        %Y(:,obs_ind) = y;
        Y(:,obs_ind) = inSig(:,k);
    end
end

%Y = inSig;

%% Inference with fixed parameters
e_param            = struct;
e_param.M0         = M0;
e_param.P0         = P0;
e_param.Y          = Y;
e_param.mind       = mindy;
e_param.dt         = dty;
e_param.f_func     = f_func;
e_param.f_param    = f_param;
e_param.q_func     = q_func;
e_param.q_param    = q_param;
e_param.param_func = @tf_set_params;
e_param.h_func     = h_func;
e_param.h_param    = h_param;
e_param.R          = R;
e_param.t_min      = t_min;
e_param.int_method = @int_cubature;
e_param.ltr_ind    = ltr_ind;
e_param.start_ind  = 1;
e_param.ode_func   = @ode45;
e_param.isteps     = 1;
e_param.e_prior = theta_prior;


e_func_gf = @(w) E_sqrt_gf_count(w,e_param);

% Transformed true parameters
theta_true = theta;
ltr_theta = theta;
ltr_theta(ltr_ind) = log(ltr_theta(ltr_ind));

%% MAP Optimization of parameter (without gradients)
% Note: the parameter posterior is usually highly non-Gaussian and
% multimodal, so optimization might get stuck to a local mode and even
% fail numerically, depending on simulation and model settings. The

do_optim = 1;

if do_optim
    opt=optimset('GradObj','off');
    opt=optimset(opt,'TolX', 1e-3);
    opt=optimset(opt,'LargeScale', 'off');
    opt=optimset(opt,'Display', 'iter');
    %%%%%
    opt.MaxFunEvals = MaxFunEvals;
    %%%%%
    
    global count;
    count = 0;
    warning('off','all'); % too many warnings really slows it down
    w_opt = fminunc(e_func_gf, w0, opt);
    warning('on','all');
    
    [ltr_theta w_opt]; % this?
else
    w_opt = w0;
end

theta_opt = w_opt;
theta_opt(ltr_ind) = exp(theta_opt(ltr_ind));
[theta_true theta_opt] % or this?

results.theta_opt = theta_opt;
results.Y = Y;


% Interpolation to a finer grid for prettier visualization
isteps = 1;
e_param2 = e_param;
e_param2.dt = dt;
e_param2.mind = mind;
e_param2.isteps = isteps;
tgrid2 = linspace(t_min,t_max,steps*isteps);

es_func_gf = @(w) ES_sqrt_gf(w,e_param2);

% Estimates with guessed parameters
[E1,MM1,PP1,MS1,PS1] = feval(es_func_gf,ltr_theta);

% Estimates with optimized parameters
[E2,MM2,PP2,MS2,PS2] = feval(es_func_gf,w_opt);

%% Plotting estimates of output signals
color1 = [0 0 1];
color2 = [1 0 0];
color3 = [0 1 0];
xx = tgrid'+dt;
xx2 = tgrid2'+dt./isteps;

nisteps = steps*isteps;
H = h_func;
% the predicted output means:
filtOutGuess = H*MM1;
smoothOutGuess = H*MS1;
filtOutOpt = H*MM2;
smoothOutOpt = H*MS2;
% the predicted output variance:
filtOutGuessVar = zeros(d,nisteps);
smoothOutGuessVar = zeros(d,nisteps);
filtOutOptVar = zeros(d,nisteps);
smoothOutOptVar = zeros(d,nisteps);

for i1 = 1:d
    for i = 1:size(PP1,3)
        filtOutGuessVar(i1,i) = H(i1,:)*PP1(:,:,i)*H(i1,:)';
        smoothOutGuessVar(i1,i) = H(i1,:)*PS1(:,:,i)*H(i1,:)';
        filtOutOptVar(i1,i) = H(i1,:)*PP2(:,:,i)*H(i1,:)';
        smoothOutOptVar(i1,i) = H(i1,:)*PS2(:,:,i)*H(i1,:)';
    end
end
figure(5); clf;
for i = 1:d
    subplot(d,2,2*i-1);
    fill([xx2' fliplr(xx2')], [(filtOutOpt(i,:)+1.96*sqrt(abs(filtOutOptVar(i,:)))) ...
        fliplr((filtOutOpt(i,:)-1.96*sqrt(abs(filtOutOptVar(i,:)))))], color1, 'edgecolor',color1);
    hold on
    fill([xx2' fliplr(xx2')], [(filtOutGuess(i,:)+1.96*sqrt(abs(filtOutGuessVar(i,:)))) ...
        fliplr((filtOutGuess(i,:)-1.96*sqrt(abs(filtOutGuessVar(i,:)))))], color2, 'edgecolor',color2); hold on;
    alpha(0.2);
    
    plot(xx2,filtOutOpt(i,:),'color',color1,'LineWidth',2);
    plot(xx2,filtOutGuess(i,:),'color',color2,'LineWidth',2);
    plot(xx,inSig(i,:),'k-','LineWidth',1);
    title(sprintf('Filtered estimate of output %d',i))
    hold off;
    xlim([t_min t_max])
    yl = ylim;
    
    subplot(d,2,2*i);
    fill([xx2' fliplr(xx2')], [(smoothOutOpt(i,:)+1.96*sqrt(abs(smoothOutOptVar(i,:)))) ...
        fliplr((smoothOutOpt(i,:)-1.96*sqrt(abs(smoothOutOptVar(i,:)))))], color1, 'edgecolor',color1);
    hold on
    fill([xx2' fliplr(xx2')], [(smoothOutGuess(i,:)+1.96*sqrt(abs(smoothOutGuessVar(i,:)))) ...
        fliplr((smoothOutGuess(i,:)-1.96*sqrt(abs(smoothOutGuessVar(i,:)))))], color2, 'edgecolor',color2); hold on;
    alpha(0.2);
    
    h1=plot(xx2,smoothOutOpt(i,:),'color',color1,'LineWidth',2);
    h2=plot(xx2,smoothOutGuess(i,:),'color',color2,'LineWidth',2);
    plot(xx,inSig(i,:),'k-','LineWidth',1);
    title(sprintf('Smoothed estimate of output %d',i))
    if i == d
        legend([h1;h2],'Optimised parameters','Guessed parameters')
    end
    hold off;
    xlim([t_min t_max])
    ylim(yl);
    
end

results.predMean = smoothOutOpt;
results.predVar = smoothOutOptVar;

% my plots
figure(7);clf;
subplot(511)
plot(modeMag')
title('Actual data')
subplot(512)
plot(Y(:,:)')
title('Scaled & downsampled outputs')
grid on
subplot(513)
plot(results.predMean(:,:)')
title('Predictive mean - Nonlinear Latent Force Model')
grid on

    
%% Plotting of latent forces


Hlf = [zeros(1,d) f_param.H];
% the predicted latent function means:
filtLFGuess = Hlf*MM1;
smoothLFGuess = Hlf*MS1;
filtLFOpt = Hlf*MM2;
smoothLFOpt = Hlf*MS2;
% the predicted latent function variance:
filtLFGuessVar = zeros(1,steps*isteps);
smoothLFGuessVar = zeros(1,steps*isteps);
filtLFOptVar = zeros(1,steps*isteps);
smoothLFOptVar = zeros(1,steps*isteps);
    
for i = 1:size(PP1,3)
    filtLFGuessVar(i) = Hlf*PP1(:,:,i)*Hlf';
    smoothLFGuessVar(i) = Hlf*PS1(:,:,i)*Hlf';
    filtLFOptVar(i) = Hlf*PP2(:,:,i)*Hlf';
    smoothLFOptVar(i) = Hlf*PS2(:,:,i)*Hlf';
end

figure(6); clf;
subplot(2,1,1);
fill([xx2' fliplr(xx2')], [(filtLFOpt+1.96*sqrt(abs(filtLFOptVar))) ...
    fliplr((filtLFOpt-1.96*sqrt(abs(filtLFOptVar))))], color1, 'edgecolor',color1);
hold on
fill([xx2' fliplr(xx2')], [(filtLFGuess+1.96*sqrt(abs(filtLFGuessVar))) ...
    fliplr((filtLFGuess-1.96*sqrt(abs(filtLFGuessVar))))], color2, 'edgecolor',color2); hold on;
alpha(0.2);

plot(xx2,filtLFOpt,'color',color1,'LineWidth',2);
plot(xx2,filtLFGuess,'color',color2,'LineWidth',2);
title('Filtered estimate of force')
hold off;
xlim([t_min t_max])
yl = ylim;

subplot(2,1,2);
fill([xx2' fliplr(xx2')], [(smoothLFOpt+1.96*sqrt(abs(smoothLFOptVar))) ...
    fliplr((smoothLFOpt-1.96*sqrt(abs(smoothLFOptVar))))], color1, 'edgecolor',color1);
hold on
fill([xx2' fliplr(xx2')], [(smoothLFGuess+1.96*sqrt(abs(smoothLFGuessVar))) ...
    fliplr((smoothLFGuess-1.96*sqrt(abs(smoothLFGuessVar))))], color2, 'edgecolor',color2); hold on;
alpha(0.2);

h1=plot(xx2,smoothLFOpt,'color',color1,'LineWidth',2);
h2=plot(xx2,smoothLFGuess,'color',color2,'LineWidth',2);
title('Smoothed estimate of force')
legend([h1;h2],'Optimised parameters','Guessed parameters')
hold off;
xlim([t_min t_max])
ylim(yl);

results.u = smoothLFOpt;
results.uVar = smoothLFOptVar;
results.g = log(1+exp(smoothLFOpt(1,:)));

results.allMag = allMag;
results.allFreq = allFreq;
results.medFreq = medFreq;

results.modeMag = modeMag;
results.modeFreq = modeFreq;
results.modeMedFreq = modeMedFreq;

results.residualMag = residualMag;
results.residualFreq = residualFreq;
results.residualMedFreq = residualMedFreq;

results.dt = dt;
results.dty = dty;
results.sfactor = sfactor;
results.maxOfMode = maxOfMode;
results.gamma = gamma;
results.polyOrder = polyOrder;


%% run the state space model

A = results.theta_opt(3:numModes+2);
A_ = A;
B = results.theta_opt(numModes+3:(2*numModes)+2);
B_ = B;
D = eye(numModes,numModes);
%D(eye(numModes)==1) = -((1/2)+results.theta_opt((2*numModes)+3:(3*numModes)+2).^2);
D(eye(numModes)==1) = -results.theta_opt((2*numModes)+3:(3*numModes)+2);
D_ = -diag(D);
S = results.theta_opt((3*numModes)+3:(4*numModes)+2);
    normFactor = max(results.g);
    g = results.g ./ normFactor;
    %g(end-15:end) = 0;
S_ = S;
S = normFactor .* S;

% run the state space model
    T = length(g); % this can be whatever really, but we'll need a latent force of the same length

    y = zeros(numModes,T);    % Preallocate output signal for t=1:T
    
    u = g; % the latent force, just use this for now, since it's the right length

    % Perform the system simulation:
    x = A;                           % Set initial state
    for t=1:T                           % Iterate through time
        y(:,t) = x;% .* adjustFactor';                 % Output for time t
        Su = S.*u(t);            % Latent force impact
        xdot = D*real(x.^(gamma)) + Su + B;  % calculate x'(t)
        x = x + (xdot * dt);        % calculate x(t+1) = forward one time step
    end

% get the right scaling variables to obtain overall magnitude close to original data 
scaleFactor = results.maxOfMode ./ max(y');

results.scaleFactor = scaleFactor;

yScaled = bsxfun(@times,y,scaleFactor');

% my plots
figure(7)
s5=subplot(515); cla(s5)
plot(smoothLFOpt(1,:))
hold on
plot(log(1+exp(smoothLFOpt(1,:))))
title('Smoothed Latent Force')
s4=subplot(514); cla(s4)
plot(yScaled')
title('State space model outputs (re-scaled)')

%% save
if saveResults == 1
    save(matNameMag,'results');
end
