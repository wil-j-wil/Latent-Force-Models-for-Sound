% this code is messy but it produces the exact graphs from my DAFx-2017
% paper, including formatting


%% Figure 1 - modelChoice
% illustrate the impact of different choices of gamma in graph form

alpha_=0.162; % sensitivity
beta_=88; % damping

T=45;
hopSize=441;
dt=1/hopSize;
dt_=0.1;

t=linspace(1,T,T);
xE=alpha_*exp(-beta_*t/hopSize);


 
 % get the real clarinet data
 results = load('resultsClarinet');
 results = results.results;
 clarMag = results.spearMag;

x=xE(1);
x_3quarters=xE(1);
x_half=xE(1);
x_quarter=xE(1);
y=zeros(T,4);
D=beta_;
D_=beta_*dt/dt_;
%D_=2.4;
D_3quarters=1.01;
D_half=0.537;
D_quarter=0.275;
for t_=1:T                           % Iterate through time
        y(t_,:) = [max(x,0) max(x_3quarters,0) max(x_half,0) max(x_quarter,0)]; % Output for time t
        xdot = -D_*real(x^(1));  % calculate x'(t)
        xdot_3quarters = -D_3quarters*real(x_3quarters^(3/4));  % calculate x'(t)
        xdot_half = -D_half*real(x_half^(1/2));  % calculate x'(t)
        xdot_quarter = -D_quarter*real(x_quarter^(1/4));  % calculate x'(t)
        x = x + (xdot * dt_);
        x_3quarters = x_3quarters + (xdot_3quarters * dt_);        % calculate x(t+1) = forward one time step
        x_half = x_half + (xdot_half * dt_);        % calculate x(t+1) = forward one time step
        x_quarter = x_quarter + (xdot_quarter * dt_);        % calculate x(t+1) = forward one time step
end
fig1=figure(91);clf
set(fig1, 'Position', [1, 650, 1*300, 1*170]);
ax1=axes('Parent',fig1,'fontsize',12,'TickLabelInterpreter', 'latex');
%subplot(212);cla
startingPoint=27;
plot(ax1,([0:T-startingPoint]+startingPoint)/100,clarMag(2,startingPoint:T),'k','LineWidth',2) % this is from the flute
hold on
plot(ax1,(t(1:T-startingPoint+1)+startingPoint-1)/100,y(1:T-startingPoint+1,1),'k-','LineWidth',0.8) % gamma = 1
plot(ax1,(t(1:T-startingPoint+1)+startingPoint-1)/100,y(1:T-startingPoint+1,2),'k--','LineWidth',0.8,'MarkerSize',2.2) % gamma = 3/4
plot(ax1,(t(1:T-startingPoint+1)+startingPoint-1)/100,y(1:T-startingPoint+1,3),'ko','LineWidth',0.8,'MarkerSize',2.2) % gamma = 1/2
%plot(ax1,(t(1:T-startingPoint+1)+startingPoint-1)/100,y(1:T-startingPoint+1,4),'ko','LineWidth',1,'MarkerSize',3) % gamma = 1/4
axis([0.268 0.435, 0 0.142])

%plot(4*modeMag(8,1:T))
l=legend(ax1,'Clarinet recording (decay)','$\gamma=1$','$\gamma=3/4$','$\gamma=1/2$');%,'$\gamma=1/4$');
set(l,'interpreter', 'Latex')
set(gca,'TickLabelInterpreter', 'latex')
%set(leg,'FontSize',11);
xlabel('Time (seconds)', 'interpreter', 'Latex');%,'FontSize',11)
ylabel('Amplitude', 'interpreter', 'Latex');%,'FontSize',11)


%% Figure 2 - modePicking / modePickingArrows
[allMag, allFreq] = readSpearFile('metal-impact.txt'); % sinusoidal analysis data from Spear
 
 % normalise so that maximum value is 0.3
 allMagMax = max(allMag');
 maxMax = max(allMagMax);
 allMag = 0.3 * allMag / maxMax;
 
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
 
 % fit a polynomial to our frequency data
 warning('off','all');
  sinPoly2 = polyfit(binInd,mag2db(binSum)',12); % poly of order 10 or 12 works well
 warning('on','all');
 sinPolyDraw2 = polyval(sinPoly2, binInd); % this will become our reference
 
 % "filter" the data
 binFilt = mag2db(binSum)-sinPolyDraw2';
 
 % plot for comparison
 fig2=figure(92); clf
 set(fig2, 'Position', [302, 530, 300, 340]);
 ax2=axes('Parent',fig2,'fontsize',12,'TickLabelInterpreter', 'latex');
 subplot(211)
 plot(binInd,mag2db(binSum),'k-')
 axis([0 25000, -80 30])
 hold on
 polyPlot=plot(binInd,sinPolyDraw2,'r-','LineWidth',1);
 ylabel('Power/Freq. (dB/Hz)', 'interpreter', 'Latex')
 title('a) Frequency spectrum of a metal impact sound', 'interpreter', 'Latex')
 set(gca,'Ytick',-75:25:25,'TickLabelInterpreter', 'latex')
 l=legend(polyPlot,'Polynomial fit');
 set(l,'interpreter', 'Latex')
 set(gca,'XTickLabel',[]);
 subplot(212)
 plot(binInd,binFilt-30,'k')
 axis([0 25000, -80 30])
 xlabel('Frequency (Hz)', 'interpreter', 'Latex')
 ylabel('Power/Freq. (dB/Hz)', 'interpreter', 'Latex')
 title('b) The filtered / flattened spectrum', 'interpreter', 'Latex')
 set(gca,'Ytick',-75:25:25,'TickLabelInterpreter', 'latex')
 %annotation('textarrow',[0.5,0.4],[0.4,0.35],'String','mode')
 annotation('arrow',[0.305,0.19],[0.415,0.38],'Color','b','HeadLength',5,'HeadWidth',4)
 annotation('arrow',[0.33,0.285],[0.41,0.38],'Color','b','HeadLength',5,'HeadWidth',4)
 annotation('arrow',[0.39,0.44],[0.41,0.37],'Color','b','HeadLength',5,'HeadWidth',4)
 annotation('textbox',[0.3 0.35 0.1 0.1],'String','modes','FitBoxToText','on', 'interpreter', 'Latex', 'EdgeColor', 'none','Color','b')
 
 
 
%% Figure 3 - resultsClarinet

resultsC = load('resultsClarinet'); % clarinet results
uVar = resultsC.results.uVar;
resultsC = resultsC.results;

% unwrap the parameters from theta vector
dt = 0.1; % should be standard across models

numModes=6;
A = resultsC.theta_opt(3:numModes+2);
B = resultsC.theta_opt(numModes+3:(2*numModes)+2);
D = eye(numModes,numModes);
D(eye(numModes)==1) = -resultsC.theta_opt((2*numModes)+3:(3*numModes)+2);
S = resultsC.theta_opt((3*numModes)+3:(4*numModes)+2);

maxMean = max(resultsC.predMean');
adjustFactor = resultsC.maxOfMode ./ maxMean;

T = length(resultsC.g); % this can be whatever really, but we'll need a latent force of the same length

y = zeros(numModes,T);    % Preallocate output signal for t=1:T

u = resultsC.g; % the latent force, just use this for now, since it's the right length

% Perform the system simulation:
x = A;                           % Set initial state
for t=1:T                           % Iterate through time
    y(:,t) = max(x, 0); % Output for time t
    Su = S.*u(t);            % Latent force impact
    xdot = D*real(x.^(1/2)) + Su + B;  % calculate x'(t)
    x = x + (xdot * dt);        % calculate x(t+1) = forward one time step
end

reScale = resultsC.maxOfMode ./ max(y');

freqC = resultsC.spearFreq(:,1:end-9);
freqC(find(freqC==0))=NaN;


fig4=figure(94);clf
set(fig4, 'Position', [603, 470, 680, 340]);
ax4=axes('Parent',fig4,'fontsize',12,'TickLabelInterpreter', 'latex');
subplot(221)
ghost1=plot(NaN,'k-');
hold on
ghost2=plot(NaN,'k--');
plot([1:43]/100,resultsC.spearMag(:,1:end-9)','LineWidth',1)
ax=gca;
ax.ColorOrderIndex=1;
plot([1:43]/100,(bsxfun(@times,resultsC.predMean(:,1:end-9),adjustFactor'))','--','LineWidth',1)
ghost1=plot(NaN,'k-');
ghost2=plot(NaN,'k--');
axis([0 0.44, 0 0.3])
set(gca,'TickLabelInterpreter', 'latex')
title('a) Predictive mean of the outputs for a clarinet note', 'interpreter', 'Latex')
ylabel('Amplitude', 'interpreter', 'Latex')
l=legend('Actual data','Pred. mean');
set(l,'interpreter', 'Latex')
set(gca,'XTickLabel',[]);
%subplot(412)
%plot((bsxfun(@times,resultsC.predMean(:,1:end-9),adjustFactor'))','LineWidth',1)
%axis([0 44, 0 0.3])
%title('Predictive mean of the outputs', 'interpreter', 'Latex')
%xlabel('Frame Number', 'interpreter', 'Latex')
%ylabel('Amplitude', 'interpreter', 'Latex')
subplot(224)
p1=fill([[1:43]/100 fliplr([1:43]/100)],[(resultsC.u(1:end-9)+1.96*sqrt(uVar(1:end-9))) ...
    fliplr((resultsC.u(1:end-9)-1.96*sqrt(uVar(1:end-9))))], [0.85 0.85 0.85], 'edgecolor',[0.85 0.85 0.85]);
hold on
p2=plot([1:43]/100,resultsC.u(1:end-9),'k--','LineWidth',1);
p3=plot([1:43]/100,resultsC.g(1:end-9),'k-','LineWidth',1);
axis([0 0.44, -5 2])
set(gca,'TickLabelInterpreter', 'latex')
title('d) Predictive mean of the input', 'interpreter', 'Latex')
l=legend([p2 p1 p3],'Latent force mean, $\bar{u}$','Uncertainty','Softplus mean, $g(\bar{u})$','Location','southwest');
set(l,'interpreter', 'Latex')
xlabel('Time (seconds)', 'interpreter', 'Latex')
subplot(222)
plot([1:43]/100,bsxfun(@times,y(:,1:end-9)',reScale),'LineWidth',1)
axis([0 0.44, 0 0.3])
set(gca,'TickLabelInterpreter', 'latex')
title('b) Resynthesis with the state space model', 'interpreter', 'Latex')
ylabel('Amplitude', 'interpreter', 'Latex')
set(gca,'XTickLabel',[]);

subplot(223)
plot([1:43]/100,freqC','LineWidth',1)
axis([0 0.44, 0 2000])
set(gca,'TickLabelInterpreter', 'latex')
title('c) Frequency of the modes', 'interpreter', 'Latex')
xlabel('Time (seconds)', 'interpreter', 'Latex')
ylabel('Frequency (Hz)', 'interpreter', 'Latex')



%% Figure 4 - resultsMetal

resultsM = load('resultsMetalImpact');
uVar = resultsM.results.uVar;
resultsM = resultsM.results;

% unwrap the parameters from theta vector
dt = 0.1; % should be standard across models

numModes=8;
A = resultsM.theta_opt(3:numModes+2);
B = resultsM.theta_opt(numModes+3:(2*numModes)+2);
D = eye(numModes,numModes);
D(eye(numModes)==1) = -resultsM.theta_opt((2*numModes)+3:(3*numModes)+2);
S = resultsM.theta_opt((3*numModes)+3:(4*numModes)+2);

maxMean = max(resultsM.predMean');
adjustFactor = resultsM.maxOfMode ./ maxMean;

T = length(resultsM.g); % this can be whatever really, but we'll need a latent force of the same length

y = zeros(numModes,T);    % Preallocate output signal for t=1:T

u = resultsM.g; % the latent force, just use this for now, since it's the right length

% Perform the system simulation:
x = A;                           % Set initial state
for t=1:T                           % Iterate through time
    y(:,t) = max(x, 0); % Output for time t
    Su = S.*u(t);            % Latent force impact
    xdot = D*real(x.^(3/4)) + Su + B;  % calculate x'(t)
    x = x + (xdot * dt);        % calculate x(t+1) = forward one time step
end

reScale = resultsM.maxOfMode ./ max(y');


fig5=figure(95);clf
set(fig5, 'Position', [0, 53, 680, 340]);
ax5=axes('Parent',fig5,'fontsize',12,'TickLabelInterpreter', 'latex');
subplot(221)
ghost1=plot(NaN,'k-');
hold on
ghost2=plot(NaN,'k--');
plot([1:92]/100,resultsM.modeMag','LineWidth',1)
ax=gca;
ax.ColorOrderIndex=1;
%plot((bsxfun(@times,resultsM.predMean,adjustFactor'))','--','LineWidth',1)
ghost1=plot(NaN,'k-');
%ghost2=plot(NaN,'k--');
axis([0 0.93, -0.035 0.35])
set(gca,'TickLabelInterpreter', 'latex')
title('a) Amplitude of 8 modes of a metal impact sound', 'interpreter', 'Latex')
ylabel('Amplitude', 'interpreter', 'Latex')
set(gca,'XTickLabel',[]);
%l=legend('Actual data','Pred. mean');
%set(l,'interpreter', 'Latex')
subplot(224)
p1=fill([[1:T]/100 fliplr([1:T]/100)],[(resultsM.u+1.96*sqrt(uVar)) ...
    fliplr((resultsM.u-1.96*sqrt(uVar)))], [0.85 0.85 0.85], 'edgecolor',[0.85 0.85 0.85]);
hold on
p2=plot([1:92]/100,resultsM.u,'k--','LineWidth',1);
p3=plot([1:92]/100,resultsM.g,'k-','LineWidth',1);
axis([0 0.93, -5 2])
set(gca,'TickLabelInterpreter', 'latex')
title('d) Predictive mean of the input', 'interpreter', 'Latex')
l=legend([p2 p1 p3],'Latent force mean, $\bar{u}$','Uncertainty','Softplus mean, $g(\bar{u})$','Location','east');
set(l,'interpreter', 'Latex')
xlabel('Time (seconds)', 'interpreter', 'Latex')
subplot(222)
plot([1:92]/100,bsxfun(@times,y',reScale),'LineWidth',1)
axis([0 0.93, -0.035 0.35])
set(gca,'TickLabelInterpreter', 'latex')
title('b) Resynthesis with the state space model', 'interpreter', 'Latex')
set(gca,'XTickLabel',[]);


% Figure 4c - pcaMetal

results = load('resultsMetalImpact');
results = results.results;

Xp = results.modeMag';

[w, pc, ev] = pca(Xp);

mu = mean(Xp);
xhat = bsxfun(@minus,Xp,mu); % subtract the mean
norm(pc * w' - xhat);

Xnew1 = pc(:,1) * w(:,1)';
Xnew1 = bsxfun(@plus,mu,Xnew1); % add the mean back in

Xnew2 = pc(:,1:2) * w(:,1:2)';
Xnew2 = bsxfun(@plus,mu,Xnew2); % add the mean back in

PCA_error1 = mean(mean((Xp-Xnew1).^2))^(1/2); % RMSE
PCA_error2 = mean(mean((Xp-Xnew2).^2))^(1/2); % RMSE

subplot(223)
plot([1:92]/100,Xnew1,'LineWidth',1)
hold on
axis([0 0.93, -0.035 0.35])
set(gca,'TickLabelInterpreter', 'latex')
title('c) Resynthesis with 1 principal component', 'interpreter', 'Latex')
xlabel('Time (seconds)', 'interpreter', 'Latex')
ylabel('Amplitude', 'interpreter', 'Latex')



%% Figure 5 - morph

interp = 0.5;
force = 0;

% load data and plot
numModels = 2;
model = cell(1,numModels);
results1 = load('resultsOboe');
model{1} = results1.results;
model{1}.freq = model{1}.spearFreq(:,10); % manually set to 10 rather than 1

results2 = load('resultsClarinet');
model{2} = results2.results;
model{2}.freq = model{2}.spearFreq(:,1);


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

% unwrap the parameters from theta vector
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
    S{j} = normFactor .* S{j};
end


% run the state space model
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
        xdot = D{m}*real(x.^(1/2)) + Su + B{m};  % calculate x'(t)
        x = x + (xdot * dt);        % calculate x(t+1) = forward one time step
    end
end

y1=y{1};
y2=y{2};

% morphing / interpolation?
%interp = 0.25; % [0,1], 0 = model1, 1 = model2
i0 = 1-interp;
i1 = interp;

m = force+1; % which latent force to use

lf = dlmread('drawnLatentForce',',');

%u = model{m}.g; % the latent force, just use this for now, since it's the right length
u = lf;
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
        Su = (i0.*S{1}+i1.*S{2}).*u(t);                         % Latent force impact)
        xdot = (i0.*D{1}+i1.*D{2})*real(x.^(1/2)) + Su + i0.*B{1}+i1.*B{2};  % calculate x'(t)
        x = x + (xdot * dt);                                    % calculate x(t+1) = forward one time step
        yFreq(:,t) = i0.*model{1}.freq + i1.*model{2}.freq;
end

freq1=model{1}.spearFreq;
freq1(find(freq1==0))=NaN;
freq2=model{2}.spearFreq;
freq2(find(freq2==0))=NaN;


fig7=figure(97); clf
set(fig7, 'Position', [682, 53, 680, 340]);
ax7=axes('Parent',fig7,'fontsize',12,'TickLabelInterpreter', 'latex');
subplot(331)
plot([1:120]/100,y1','LineWidth',1)
axis([0 1.20, 0 0.35])
set(gca,'Xtick',0:0.40:1.20,'TickLabelInterpreter', 'latex')
%title('Amplitude of modes of an oboe', 'interpreter', 'Latex')
title('a) Oboe', 'interpreter', 'Latex')
ylabel('Amplitude', 'interpreter', 'Latex')
set(gca,'XTickLabel',[]);

subplot(334)
plot([1:120]/100,model{1}.g,'k-','LineWidth',1)
axis([0 1.20, 0 1])
set(gca,'Xtick',0:.40:1.20,'TickLabelInterpreter', 'latex')
%title('Oboe excitation function', 'interpreter', 'Latex')
ylabel('Excitation function', 'interpreter', 'Latex')
set(gca,'XTickLabel',[]);

subplot(333)
plot([1:52]/100,y2','LineWidth',1)
axis([0 0.50, 0 0.35])
set(gca,'Xtick',0:.20:.40,'TickLabelInterpreter', 'latex')
%title('Amplitude of modes of a clarinet', 'interpreter', 'Latex')
title('c) Clarinet', 'interpreter', 'Latex')
%ylabel('Amplitude', 'interpreter', 'Latex')
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);

subplot(336)
plot([1:52]/100,model{2}.g,'k-','LineWidth',1)
axis([0 .50, 0 1])
set(gca,'Xtick',0:.20:.40,'TickLabelInterpreter', 'latex')
%title('Clarinet excitation function', 'interpreter', 'Latex')
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);

subplot(332)
plot([1:142]/100,y','LineWidth',1)
axis([0 1.00, 0 0.35])
set(gca,'Xtick',0:.40:1.60,'TickLabelInterpreter', 'latex')
%title('Amplitude of morphed modes', 'interpreter', 'Latex')
title('b) Morph', 'interpreter', 'Latex')
%ylabel('Amplitude', 'interpreter', 'Latex')
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);

subplot(335)
plot([1:142]/100,u,'k-','LineWidth',1)
axis([0 1.00, 0 1])
set(gca,'Xtick',0:.40:1.60,'TickLabelInterpreter', 'latex')
%title('Excitation function used for morph', 'interpreter', 'Latex')
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);

subplot(337)
plot([1:120]/100,freq1','LineWidth',1)
axis([0 1.20, 0 2000])
set(gca,'Xtick',0:.40:1.20,'TickLabelInterpreter', 'latex')
%title('Frequency of modes of an oboe', 'interpreter', 'Latex')
xlabel('Time (seconds)', 'interpreter', 'Latex')
ylabel('Frequency (Hz)', 'interpreter', 'Latex')

subplot(339)
plot([1:52]/100,freq2(1,:)','LineWidth',1)
hold on
ax=gca;
ax.ColorOrderIndex=3;
plot([1:52]/100,freq2(2:end,:)','LineWidth',1)
axis([0 .50, 0 2000])
set(gca,'Xtick',0:.20:.40,'TickLabelInterpreter', 'latex')
%title('Frequency of modes of a clarinet', 'interpreter', 'Latex')
xlabel('Time (seconds)', 'interpreter', 'Latex')
set(gca,'YTickLabel',[]);

subplot(338)
plot([1:142]/100,yFreq','LineWidth',1)
axis([0 1.00, 0 2000])
set(gca,'Xtick',0:.40:1.20,'TickLabelInterpreter', 'latex')
%title('Frequency of morphed modes', 'interpreter', 'Latex')
xlabel('Time (seconds)', 'interpreter', 'Latex')
set(gca,'YTickLabel',[]);



%% table 1 - RMSE

% metal
results = load('resultsMetalImpact');
results = results.results;

% LFM
% unwrap the parameters from theta vector
dt = 0.1; % should be standard across models

numModes=8;
A = results.theta_opt(3:numModes+2);
B = results.theta_opt(numModes+3:(2*numModes)+2);
D = eye(numModes,numModes);
D(eye(numModes)==1) = -results.theta_opt((2*numModes)+3:(3*numModes)+2);
S = results.theta_opt((3*numModes)+3:(4*numModes)+2);

maxMean = max(results.predMean');
adjustFactor = results.maxOfMode ./ maxMean;

T = length(results.g); % this can be whatever really, but we'll need a latent force of the same length

y = zeros(numModes,T);    % Preallocate output signal for t=1:T

u = results.g; % the latent force, just use this for now, since it's the right length

% Perform the system simulation:
x = A;                           % Set initial state
for t=1:T                           % Iterate through time
    y(:,t) = max(x, 0); % Output for time t
    Su = S.*u(t);            % Latent force impact
    xdot = D*real(x.^(3/4)) + Su + B;  % calculate x'(t)
    x = x + (xdot * dt);        % calculate x(t+1) = forward one time step
end

reScale = results.maxOfMode ./ max(y');
LFMy = bsxfun(@times,y',reScale);

% PCA
Xp = results.modeMag';

[w, pc, ev] = pca(Xp);

mu = mean(Xp);
xhat = bsxfun(@minus,Xp,mu); % subtract the mean
norm(pc * w' - xhat);

Xnew1 = pc(:,1) * w(:,1)';
Xnew1 = bsxfun(@plus,mu,Xnew1); % add the mean back in

Xnew2 = pc(:,1:2) * w(:,1:2)';
Xnew2 = bsxfun(@plus,mu,Xnew2); % add the mean back in

muL = mean(LFMy);
xhatL = bsxfun(@minus,LFMy,muL); % subtract the mean
muP1 = mean(Xnew1);
xhatP1 = bsxfun(@minus,Xnew1,muP1); % subtract the mean
muP2 = mean(Xnew2);
xhatP2 = bsxfun(@minus,Xnew2,muP2); % subtract the mean

% normalise
normFactor=0.3./max(abs(xhat));
x_=bsxfun(@times,xhat,normFactor);
normFactor=0.3./max(abs(xhatL));
x_L=bsxfun(@times,xhatL,normFactor);
normFactor=0.3./max(abs(xhatP1));
x_P1=bsxfun(@times,xhatP1,normFactor);
normFactor=0.3./max(abs(xhatP2));
x_P2=bsxfun(@times,xhatP2,normFactor);

LFM_errorM  = mean(mean((x_-x_L).^2))^(1/2); % RMSE
PCA_errorM1 = mean(mean((x_-x_P1).^2))^(1/2); % RMSE
PCA_errorM2 = mean(mean((x_-x_P2).^2))^(1/2); % RMSE


% wood
results = load('resultsWoodenImpact_2'); % this was a seperate run with gamma = 1
results = results.results;

% LFM
% unwrap the parameters from theta vector
dt = 0.1; % should be standard across models

numModes=8;
A = results.theta_opt(3:numModes+2);
B = results.theta_opt(numModes+3:(2*numModes)+2);
D = eye(numModes,numModes);
D(eye(numModes)==1) = -results.theta_opt((2*numModes)+3:(3*numModes)+2);
S = results.theta_opt((3*numModes)+3:(4*numModes)+2);

maxMean = max(results.predMean');
adjustFactor = results.maxOfMode ./ maxMean;

T = length(results.g); % this can be whatever really, but we'll need a latent force of the same length

y = zeros(numModes,T);    % Preallocate output signal for t=1:T

u = results.g; % the latent force, just use this for now, since it's the right length

% Perform the system simulation:
x = A;                           % Set initial state
for t=1:T                           % Iterate through time
    y(:,t) = max(x, 0); % Output for time t
    Su = S.*u(t);            % Latent force impact
    %xdot = D*real(x.^(3/4)) + Su + B;  % calculate x'(t)
    xdot = D*x + Su + B;  % calculate x'(t)
    x = x + (xdot * dt);        % calculate x(t+1) = forward one time step
end

reScale = results.maxOfMode ./ max(y');
LFMy = bsxfun(@times,y',reScale);

% PCA
Xp = results.modeMag';

[w, pc, ev] = pca(Xp);

mu = mean(Xp);
xhat = bsxfun(@minus,Xp,mu); % subtract the mean
norm(pc * w' - xhat);

Xnew1 = pc(:,1) * w(:,1)';
Xnew1 = bsxfun(@plus,mu,Xnew1); % add the mean back in

Xnew2 = pc(:,1:2) * w(:,1:2)';
Xnew2 = bsxfun(@plus,mu,Xnew2); % add the mean back in

muL = mean(LFMy);
xhatL = bsxfun(@minus,LFMy,muL); % subtract the mean
muP1 = mean(Xnew1);
xhatP1 = bsxfun(@minus,Xnew1,muP1); % subtract the mean
muP2 = mean(Xnew2);
xhatP2 = bsxfun(@minus,Xnew2,muP2); % subtract the mean

% normalise
normFactor=0.3./max(abs(xhat));
x_=bsxfun(@times,xhat,normFactor);
normFactor=0.3./max(abs(xhatL));
x_L=bsxfun(@times,xhatL,normFactor);
normFactor=0.3./max(abs(xhatP1));
x_P1=bsxfun(@times,xhatP1,normFactor);
normFactor=0.3./max(abs(xhatP2));
x_P2=bsxfun(@times,xhatP2,normFactor);

LFM_errorW  = mean(mean((x_-x_L).^2))^(1/2); % RMSE
PCA_errorW1 = mean(mean((x_-x_P1).^2))^(1/2); % RMSE
PCA_errorW2 = mean(mean((x_-x_P2).^2))^(1/2); % RMSE





% clarinet
results = load('resultsClarinet');
results = results.results;

% LFM
% unwrap the parameters from theta vector
dt = 0.1; % should be standard across models

numModes=6;
A = results.theta_opt(3:numModes+2);
B = results.theta_opt(numModes+3:(2*numModes)+2);
D = eye(numModes,numModes);
D(eye(numModes)==1) = -results.theta_opt((2*numModes)+3:(3*numModes)+2);
S = results.theta_opt((3*numModes)+3:(4*numModes)+2);

maxMean = max(results.predMean');
adjustFactor = results.maxOfMode ./ maxMean;

T = length(results.g)-9; % this can be whatever really, but we'll need a latent force of the same length

y = zeros(numModes,T);    % Preallocate output signal for t=1:T

u = results.g; % the latent force, just use this for now, since it's the right length

% Perform the system simulation:
x = A;                           % Set initial state
for t=1:T                           % Iterate through time
    y(:,t) = max(x, 0); % Output for time t
    Su = S.*u(t);            % Latent force impact
    xdot = D*real(x.^(1/2)) + Su + B;  % calculate x'(t)
    x = x + (xdot * dt);        % calculate x(t+1) = forward one time step
end

reScale = results.maxOfMode ./ max(y');
LFMy = bsxfun(@times,y',reScale);

% PCA
Xp = results.spearMag(:,1:end-9)';

[w, pc, ev] = pca(Xp);

mu = mean(Xp);
xhat = bsxfun(@minus,Xp,mu); % subtract the mean
norm(pc * w' - xhat);

Xnew1 = pc(:,1) * w(:,1)';
Xnew1 = bsxfun(@plus,mu,Xnew1); % add the mean back in

Xnew2 = pc(:,1:2) * w(:,1:2)';
Xnew2 = bsxfun(@plus,mu,Xnew2); % add the mean back in

muL = mean(LFMy);
xhatL = bsxfun(@minus,LFMy,muL); % subtract the mean
muP1 = mean(Xnew1);
xhatP1 = bsxfun(@minus,Xnew1,muP1); % subtract the mean
muP2 = mean(Xnew2);
xhatP2 = bsxfun(@minus,Xnew2,muP2); % subtract the mean

% normalise
normFactor=0.3./max(abs(xhat));
x_=bsxfun(@times,xhat,normFactor);
normFactor=0.3./max(abs(xhatL));
x_L=bsxfun(@times,xhatL,normFactor);
normFactor=0.3./max(abs(xhatP1));
x_P1=bsxfun(@times,xhatP1,normFactor);
normFactor=0.3./max(abs(xhatP2));
x_P2=bsxfun(@times,xhatP2,normFactor);

LFM_errorC  = mean(mean((x_-x_L).^2))^(1/2); % RMSE
PCA_errorC1 = mean(mean((x_-x_P1).^2))^(1/2); % RMSE
PCA_errorC2 = mean(mean((x_-x_P2).^2))^(1/2); % RMSE



% oboe
results = load('resultsOboe');
results = results.results;

% LFM
% unwrap the parameters from theta vector
dt = 0.1; % should be standard across models

numModes=6;
A = results.theta_opt(3:numModes+2);
B = results.theta_opt(numModes+3:(2*numModes)+2);
D = eye(numModes,numModes);
D(eye(numModes)==1) = -results.theta_opt((2*numModes)+3:(3*numModes)+2);
S = results.theta_opt((3*numModes)+3:(4*numModes)+2);

maxMean = max(results.predMean');
adjustFactor = results.maxOfMode ./ maxMean;

T = length(results.g)-9; % this can be whatever really, but we'll need a latent force of the same length

y = zeros(numModes,T);    % Preallocate output signal for t=1:T

u = results.g; % the latent force, just use this for now, since it's the right length

% Perform the system simulation:
x = A;                           % Set initial state
for t=1:T                           % Iterate through time
    y(:,t) = max(x, 0); % Output for time t
    Su = S.*u(t);            % Latent force impact
    xdot = D*real(x.^(1/2)) + Su + B;  % calculate x'(t)
    x = x + (xdot * dt);        % calculate x(t+1) = forward one time step
end

reScale = results.maxOfMode ./ max(y');
LFMy = bsxfun(@times,y',reScale);

% PCA
Xp = results.spearMag(:,1:end-9)';

[w, pc, ev] = pca(Xp);

mu = mean(Xp);
xhat = bsxfun(@minus,Xp,mu); % subtract the mean
norm(pc * w' - xhat);

Xnew1 = pc(:,1) * w(:,1)';
Xnew1 = bsxfun(@plus,mu,Xnew1); % add the mean back in

Xnew2 = pc(:,1:2) * w(:,1:2)';
Xnew2 = bsxfun(@plus,mu,Xnew2); % add the mean back in

muL = mean(LFMy);
xhatL = bsxfun(@minus,LFMy,muL); % subtract the mean
muP1 = mean(Xnew1);
xhatP1 = bsxfun(@minus,Xnew1,muP1); % subtract the mean
muP2 = mean(Xnew2);
xhatP2 = bsxfun(@minus,Xnew2,muP2); % subtract the mean

% normalise
normFactor=0.3./max(abs(xhat));
x_=bsxfun(@times,xhat,normFactor);
normFactor=0.3./max(abs(xhatL));
x_L=bsxfun(@times,xhatL,normFactor);
normFactor=0.3./max(abs(xhatP1));
x_P1=bsxfun(@times,xhatP1,normFactor);
normFactor=0.3./max(abs(xhatP2));
x_P2=bsxfun(@times,xhatP2,normFactor);

LFM_errorO  = mean(mean((x_-x_L).^2))^(1/2); % RMSE
PCA_errorO1 = mean(mean((x_-x_P1).^2))^(1/2); % RMSE
PCA_errorO2 = mean(mean((x_-x_P2).^2))^(1/2); % RMSE




% piano
results = load('resultsPiano');
results = results.results;

% LFM
% unwrap the parameters from theta vector
dt = 0.1; % should be standard across models

numModes=6;
A = results.theta_opt(3:numModes+2);
B = results.theta_opt(numModes+3:(2*numModes)+2);
D = eye(numModes,numModes);
D(eye(numModes)==1) = -results.theta_opt((2*numModes)+3:(3*numModes)+2);
S = results.theta_opt((3*numModes)+3:(4*numModes)+2);

maxMean = max(results.predMean');
adjustFactor = results.maxOfMode ./ maxMean;

T = length(results.g)-9; % this can be whatever really, but we'll need a latent force of the same length

y = zeros(numModes,T);    % Preallocate output signal for t=1:T

u = results.g; % the latent force, just use this for now, since it's the right length

% Perform the system simulation:
x = A;                           % Set initial state
for t=1:T                           % Iterate through time
    y(:,t) = max(x, 0); % Output for time t
    Su = S.*u(t);            % Latent force impact
    xdot = D*real(x.^(1/2)) + Su + B;  % calculate x'(t)
    x = x + (xdot * dt);        % calculate x(t+1) = forward one time step
end

reScale = results.maxOfMode ./ max(y');
LFMy = bsxfun(@times,y',reScale);

% PCA
Xp = results.spearMag(:,1:end-9)';

[w, pc, ev] = pca(Xp);

mu = mean(Xp);
xhat = bsxfun(@minus,Xp,mu); % subtract the mean
norm(pc * w' - xhat);

Xnew1 = pc(:,1) * w(:,1)';
Xnew1 = bsxfun(@plus,mu,Xnew1); % add the mean back in

Xnew2 = pc(:,1:2) * w(:,1:2)';
Xnew2 = bsxfun(@plus,mu,Xnew2); % add the mean back in

muL = mean(LFMy);
xhatL = bsxfun(@minus,LFMy,muL); % subtract the mean
muP1 = mean(Xnew1);
xhatP1 = bsxfun(@minus,Xnew1,muP1); % subtract the mean
muP2 = mean(Xnew2);
xhatP2 = bsxfun(@minus,Xnew2,muP2); % subtract the mean

% normalise
normFactor=0.3./max(abs(xhat));
x_=bsxfun(@times,xhat,normFactor);
normFactor=0.3./max(abs(xhatL));
x_L=bsxfun(@times,xhatL,normFactor);
normFactor=0.3./max(abs(xhatP1));
x_P1=bsxfun(@times,xhatP1,normFactor);
normFactor=0.3./max(abs(xhatP2));
x_P2=bsxfun(@times,xhatP2,normFactor);

LFM_errorP  = mean(mean((x_-x_L).^2))^(1/2); % RMSE
PCA_errorP1 = mean(mean((x_-x_P1).^2))^(1/2); % RMSE
PCA_errorP2 = mean(mean((x_-x_P2).^2))^(1/2); % RMSE

%LFM_errorP  = mean(mean((Xp-LFMy).^2))^(1/2) % RMSE
%LFM_errorP  = mean(mean((xhat-xhatL).^2))^(1/2) % RMSE
%PCA_errorP1 = mean(mean((Xp-Xnew1).^2))^(1/2) % RMSE
%PCA_errorP2 = mean(mean((Xp-Xnew2).^2))^(1/2) % RMSE

%LFM_errorP  = mean(mean((xhat-xhatL).^2))^(1/2) % RMSE
%PCA_errorP1 = mean(mean((xhat-xhatP1).^2))^(1/2) % RMSE
%PCA_errorP2 = mean(mean((xhat-xhatP2).^2))^(1/2) % RMSE

%table(LFM_errorM,PCA_errorM1,PCA_errorM2 , LFM_errorW,PCA_errorW1,PCA_errorW2 , LFM_errorC,PCA_errorC1,PCA_errorC2 ...
%    , LFM_errorO,PCA_errorO1,PCA_errorO2 , LFM_errorP,PCA_errorP1,PCA_errorP2)

LFM = [LFM_errorC; LFM_errorO; LFM_errorP; LFM_errorM; LFM_errorW];
PCA1 = [PCA_errorC1; PCA_errorO1; PCA_errorP1; PCA_errorM1; PCA_errorW1];
PCA2 = [PCA_errorC2; PCA_errorO2; PCA_errorP2; PCA_errorM2; PCA_errorW2];

RMS_data = [LFM_errorM PCA_errorM1 PCA_errorM2 ; LFM_errorW PCA_errorW1 PCA_errorW2 ; LFM_errorC PCA_errorC1 PCA_errorC2 ...
    ; LFM_errorO PCA_errorO1 PCA_errorO2 ; LFM_errorP PCA_errorP1 PCA_errorP2];
row_names = {'Clarinet','Oboe','Piano','Metal','Wood'};

RMS_error=table(LFM,PCA1,PCA2,'RowNames',row_names)