function [k2fa] = IlluminationFreqF(S1aTnoisy,OTFo)
% AIM: illumination frequency vector determination
% INPUT VARIABLES
%   S1aTnoisy: raw SIM image
%   OTFo: system OTF
% OUTPUT VARIABLE
%   k2fa: illumination frequency vector

w = size(OTFo,1);
wo = w/2;

% computing PSFe for edge tapering SIM images (heuristically determined)
PSFd = real(fftshift( ifft2(fftshift(OTFo.^10)) ));
PSFd = PSFd/max(max(PSFd));
PSFd = PSFd/sum(sum(PSFd));
h = 30;
PSFe = PSFd(wo-h+1:wo+h,wo-h+1:wo+h);

% edge tapering raw SIM image
S1aTnoisy_et = edgetaper(S1aTnoisy,PSFe);
fS1aTnoisy_et = fftshift(fft2(S1aTnoisy_et));

% OTF cut-off freq
Kotf = OTFedgeF(OTFo);

% Approx illumination frequency vector
[k2fa,~,~] = ApproxFreqDuplex(fS1aTnoisy_et,Kotf);

fS1aTnoisy = fftshift(fft2(S1aTnoisy));
% illumination frequency vector determination by optimizing
% autocorrelation of fS1aTnoisy
OPT = 1;
PhaseKai2opt0 = @(k2fa0)PhaseKai2opt(k2fa0,fS1aTnoisy,OTFo,OPT);
options = optimset('LargeScale','off','Algorithm',...
	'active-set','MaxFunEvals',500,'MaxIter',500,'Display','notify');
k2fa0 = k2fa;
[k2fa,fval] = fminsearch(PhaseKai2opt0,k2fa0,options);
% k2a = sqrt(k2fa*k2fa')
