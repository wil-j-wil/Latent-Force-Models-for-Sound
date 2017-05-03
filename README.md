# latent-force-models-for-sound
MATLAB code applying latent force modelling techniques to sinusoidal amplitude data.

This code is dependent on the LFM toolbox created by Jouni Hartikeinen and Simo Sarkka: http://becs.aalto.fi/en/research/bayes/lfm/

As input it takes export txt files created by the SPEAR sinusoidal analysis software: http://www.klingbeil.com/spear/
To apply the code to your own recordings you must download SPEAR, open your file in the software, then click file>export and save as a txt file.

Use this text file as input to the latent force model.

The main example scripts that can be edited to run on your own code are:
nl_lfmClarinet.m,
nl_lfmOboe.m,
nl_lfmMetalImpact.m,
nl_lfmWoodenImpact.m


All these scripts are the same, but with different initial settings specific to the recording they are analysing.

To perform resynthesis, see the morph scripts:
morphMusicalInstruments.m,
morphImpacts.m

running morphImpacts(0,0) synthesises the metal impact modes. morphImpacts(1,1) synthesises the wooden impact modes. morphImpacts(0.5,0) or morphImpacts(0.5,1) morphs between the two.
