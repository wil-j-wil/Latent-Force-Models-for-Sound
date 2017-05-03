function [y, yFreq] = morphMusicalInstruments(interp,force)
% slightly different to the impact one because these results were obtained
% with old LFM code and the results struct is formatted differently
if nargin < 2
    force = 0;
end

if nargin < 1
    interp = 0.5;
end

gamma = 1/2;

%% load data and plot
numModels = 2;
model = cell(1,numModels);
results1 = load('resultsOboe');
model{1} = results1.results;
model{1}.freq = model{1}.spearFreq(:,10); % manually set to 10 rather than 1

results2 = load('resultsClarinet');
model{2} = results2.results;
model{2}.freq = model{2}.spearFreq(:,1);

figure(1); clf
subplot(321)
plot(model{1}.spearMag')
title('oboe')
subplot(322)
plot(model{2}.spearMag')
title('clarinet')
subplot(323)
plot(model{1}.g)
subplot(324)
plot(model{2}.g)
subplot(325)
plot(model{1}.spearFreq')
subplot(326)
plot(model{2}.spearFreq')

% index the modes by frequency and add in zeros if a harmonic is missing
numModes = 6;
% oboe is missing none / 7th harmonic
miss1 = 7;
initPointer = 2 + miss1;
model{1}.theta_impute = model{1}.theta_opt; % the gp params
for i=1:4 % there are always 4 parameter vectors, A, B, D, S
    if i==3
        model{1}.theta_impute = [model{1}.theta_impute(1:initPointer-1); 1; model{1}.theta_impute(initPointer:end)];
    else
        model{1}.theta_impute = [model{1}.theta_impute(1:initPointer-1); 0; model{1}.theta_impute(initPointer:end)];
    end
    initPointer = initPointer + numModes + 1;
end
model{1}.freq = [model{1}.freq(1:miss1-1); 0; model{1}.freq(miss1:end)];
maxMean = max(model{1}.predMean');
model{1}.adjustFactor = model{1}.maxOfMode ./ maxMean;
model{1}.adjustFactor = [model{1}.adjustFactor(1:miss1-1)'; 0; model{1}.adjustFactor(miss1:end)'];

% clarinet is missing 2nd harmonic
miss2 = 2;
initPointer = 2 + miss2;
model{2}.theta_impute = model{2}.theta_opt; % the gp params
for i=1:4 % there are always 4 parameter vectors, A, B, D, S
    if i==3
        model{2}.theta_impute = [model{2}.theta_impute(1:initPointer-1); 1; model{2}.theta_impute(initPointer:end)];
    else
        model{2}.theta_impute = [model{2}.theta_impute(1:initPointer-1); 0; model{2}.theta_impute(initPointer:end)];
    end
    initPointer = initPointer + numModes + 1;
end
model{2}.freq = [model{2}.freq(1:miss2-1); 0; model{2}.freq(miss2:end)];
maxMean = max(model{2}.predMean');
model{2}.adjustFactor = model{2}.maxOfMode ./ maxMean;
model{2}.adjustFactor = [model{2}.adjustFactor(1:miss2-1)'; 0; model{2}.adjustFactor(miss2:end)'];

model{1}.freq(miss1) = model{2}.freq(miss1);
model{2}.freq(miss2) = model{1}.freq(miss2);
model{1}.adjustFactor(miss1) = model{2}.adjustFactor(miss1);
model{2}.adjustFactor(miss2) = model{1}.adjustFactor(miss2);

numModes = numModes + 1;

%% unwrap the parameters from theta vector
dt = 0.1; % should be standard across models

A = cell(1,numModels);
B = cell(1,numModels);
D = cell(1,numModels);
S = cell(1,numModels);
for j=1:numModels
    A{j} = model{j}.theta_impute(3:numModes+2);
    B{j} = model{j}.theta_impute(numModes+3:(2*numModes)+2);
    D{j} = eye(numModes,numModes);
    D{j}(eye(numModes)==1) = -model{j}.theta_impute((2*numModes)+3:(3*numModes)+2);
    S{j} = model{j}.theta_impute((3*numModes)+3:(4*numModes)+2);
    normFactor = max(model{j}.g);
    model{j}.g = model{j}.g ./ normFactor;
    %model{j}.g(end-15:end) = 0;
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
        y{m}(:,t) = max(x .* model{m}.adjustFactor, 0); % Output for time t
        Su = S{m}.*u(t);            % Latent force impact
        xdot = D{m}*real(x.^(gamma)) + Su + B{m};  % calculate x'(t)
        x = x + (xdot * dt);        % calculate x(t+1) = forward one time step
    end
end

figure(2); clf
subplot(231)
plot(y{1}')
title('oboe')
subplot(234)
plot(model{1}.g)

subplot(232)
plot(y{2}')
title('clarinet')
subplot(235)
plot(model{2}.g)

%% morphing / interpolation?
%interp = 0.25; % [0,1], 0 = model1, 1 = model2
i0 = 1-interp;
i1 = interp;

m = force+1; % which latent force to use

u = model{m}.g; % the latent force, just use this for now, since it's the right length
%u(40:end) = 0;
u = [u zeros(20,1)'];

T = length(u); % this can be whatever really, but we'll need a latent force of the same length

y = zeros(numModes,T);  % Preallocate output signal for t=1:T
yFreq = y; % Preallocate frequency matrix


B{1} = zeros(size(B{1}));
B{2} = zeros(size(B{2}));

ai = i0.*model{1}.adjustFactor + i1.*model{2}.adjustFactor;
% Perform the system simulation:
x = (i0.*A{1}+i1.*A{2});                                        % Set initial state
for t=1:T                                                       % Iterate through time
        y(:,t) = max(x .* ai , 0);                            % Output for time t
        Su = (i0.*S{1}+i1.*S{2}).*u(t);                         % Latent force impact
        xdot = (i0.*D{1}+i1.*D{2})*real(x.^(gamma)) + Su + i0.*B{1}+i1.*B{2};  % calculate x'(t)
        x = x + (xdot * dt);                                    % calculate x(t+1) = forward one time step
        yFreq(:,t) = i0.*model{1}.freq + i1.*model{2}.freq;
end

figure(2)
subplot(233)
plot(y')
title('morph')
subplot(236)
plot(u)

%% Play
% normalise to avoid clipping:
y = y * (0.3 / max(max(y)));

fs = 44100;
frameSize   = fs/10;
hopSize     = fs/100;
synthSound = synthtraxEdit(yFreq,y,fs,frameSize,hopSize);
sound(synthSound,fs)