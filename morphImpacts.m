function [y, yFreq] = morphImpacts(interp,force)

if nargin < 2
    force = 0;
end

if nargin < 1
    interp = 0.5;
end

gamma = 3/4;

%% load data and plot
numModels = 2;
numModes = 8;
model = cell(1,numModels);
results1 = load('resultsMetalImpact');
model{1} = results1.results;
model{1}.freq = model{1}.modeMedFreq;

results2 = load('resultsWoodenImpact');
model{2} = results2.results;
model{2}.freq = model{2}.modeMedFreq;

% make sure latent function decreases to zero at the end
model{1}.g(end-10:end) = 0;
model{2}.g(end-10:end) = 0;

% hack:
model{1}.g(40:end) = 0;

figure(101); clf
subplot(331)
plot(model{1}.modeMag')
title('metal')
subplot(332)
plot(model{2}.modeMag')
title('wood')
subplot(334)
plot(model{1}.g)
subplot(335)
plot(model{2}.g)
subplot(337)
plot(model{1}.modeFreq')
subplot(338)
plot(model{2}.modeFreq')

% calculate the adjustments

maxMean = max(model{1}.predMean');
model{1}.adjustFactor = model{1}.maxOfMode ./ maxMean;

maxMean = max(model{2}.predMean');
model{2}.adjustFactor = model{2}.maxOfMode ./ maxMean;

%% unwrap the parameters from theta vector
dt = 0.1; % should be standard across models

A = cell(1,numModels);
B = cell(1,numModels);
D = cell(1,numModels);
S = cell(1,numModels);
for j=1:numModels
    A{j} = model{j}.theta_opt(3:numModes+2);
    B{j} = model{j}.theta_opt(numModes+3:(2*numModes)+2);
    D{j} = eye(numModes,numModes);
    D{j}(eye(numModes)==1) = -model{j}.theta_opt((2*numModes)+3:(3*numModes)+2);
    S{j} = model{j}.theta_opt((3*numModes)+3:(4*numModes)+2);
    normFactor = max(model{j}.g);
    model{j}.g = model{j}.g ./ normFactor;
    S{j} = normFactor .* S{j};
end


%% run the state space model
y = cell(1,2);
for m=1:2
    T = length(model{m}.g); % this can be whatever really, but we'll need a latent force of the same length

    y{m} = zeros(numModes,T);    % Preallocate output signal for t=1:T

    u = model{m}.g; % the latent force, just use this for now, since it's the right length

    % Perform the system simulation:
    x = A{m};                           % Set initial state
    for t=1:T                           % Iterate through time
        y{m}(:,t) = max(x .* model{m}.scaleFactor', 0); % Output for time t
        Su = S{m}.*u(t);            % Latent force impact
        xdot = D{m}*real(x.^(gamma)) + Su + B{m};  % calculate x'(t)
        x = x + (xdot * dt);        % calculate x(t+1) = forward one time step
    end
end

figure(102); clf
subplot(231)
plot(y{1}')
title('metal')
subplot(234)
plot(model{1}.g)

subplot(232)
plot(y{2}')
title('wood')
subplot(235)
plot(model{2}.g)

%% morphing / interpolation?

%interp = 0.25; % [0,1], 0 = model1, 1 = model2
i0 = 1-interp;
i1 = interp;

m = force+1; % which latent force to use

u = model{m}.g; % the latent force, just use this for now, since it's the right length
%u(40:end) = 0;
%u = [u zeros(20,1)'];

T = length(u); % this can be whatever really, but we'll need a latent force of the same length

y = zeros(numModes,T);  % Preallocate output signal for t=1:T
yFreq = y; % Preallocate frequency matrix


B{1} = zeros(size(B{1}));
B{2} = zeros(size(B{2}));


ai = i0.*model{1}.scaleFactor + i1.*model{2}.scaleFactor;
% Perform the system simulation:
x = (i0.*A{1}+i1.*A{2});                                        % Set initial state
for t=1:T                                                       % Iterate through time
        y(:,t) = max(x .* ai' , 0);                            % Output for time t
        Su = (i0.*S{1}+i1.*S{2}).*u(t);                         % Latent force impact
        xdot = (i0.*D{1}+i1.*D{2})*real(x.^(gamma)) + Su + i0.*B{1}+i1.*B{2};  % calculate x'(t)
        x = x + (xdot * dt);                                    % calculate x(t+1) = forward one time step
        yFreq(:,t) = exp(i0.*log(model{1}.freq) + i1.*log(model{2}.freq)); % interpolate the (log) frequency
end

figure(102)
subplot(233)
plot(y')
title('morph')
subplot(236)
plot(u)

figure(101)
subplot(333)
plot(y')
title('morph')
subplot(336)
plot(u)
subplot(339)
plot(yFreq')

%% Play
% normalise to avoid clipping:
%y = y * (0.3 / max(max(y)));

fs = 44100;
frameSize   = fs/10;
hopSize     = fs/100;
synthSound = synthtraxEdit(yFreq,y,fs,frameSize,hopSize);
sound(synthSound,fs)