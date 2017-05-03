function [spearMag,spearFreq] = readSpearFile(fileName)
 spearFile = dlmread(fileName,'',5,0); % e.g. 'cowMoo30dB20msJoins.txt'
 % column 1 = frame start time (secs) - discard
 % column 2 = number of partials in frame - keep
 spearFile = spearFile(:,2:end);
 numFrames = size(spearFile,1); % total number of frames
 numPartials = max(spearFile(:,2))+1; % total number of partials
 
 spearMag = zeros(numPartials,numFrames);
 spearFreq = zeros(numPartials,numFrames);
 for i=1:numFrames
     framePartials = spearFile(i,1); % number of partials in this frame
     index = 2; % current position in file
     for j=1:framePartials
         spearMag(spearFile(i,index)+1,i) = spearFile(i,index+2); % magnitude
         spearFreq(spearFile(i,index)+1,i) = spearFile(i,index+1); % frequency
         index = index + 3;
     end
 end
end