function [phaseA1] = IlluminationPhaseF(S1aTnoisy,k2fa)
% AIM: illumination phase shift determination
% INPUT VARIABLES
%   S1aTnoisy: raw SIM image
%   k2fa: illumination frequency vector
% OUTPUT VARIABLE
%   phaseA1: illumination phase shift determined

PatternPhaseOpt0 = @(phaseA0)PatternPhaseOpt(phaseA0,S1aTnoisy,k2fa);
options = optimset('LargeScale','off','Algorithm',...
	'active-set','MaxFunEvals',500,'MaxIter',500,'Display','notify');
phaseA0 = 0;
[phaseA1,fval] = fminsearch(PatternPhaseOpt0,phaseA0,options);
% phaseA1*180/pi