function [fDo,fDp,fDm]...
    = PCMseparateF(S1aTnoisy,S2aTnoisy,S3aTnoisy,OTFo)
% AIM: obtaining the noisy estimates of three frequency components
% INPUT VARIABLES
%   S1aTnoisy,S2aTnoisy,S3aTnoisy: 3 raw SIM images with identical 
%                                   illumination pattern orientation 
%                                   but different phase shifts
%   OTFo: system OTF
% OUTPUT VARIABLES
%   fDo,fDp,fDm: noisy estimates of separated frequency components


w = size(S1aTnoisy,1);
wo = w/2;

% computing PSFe for edge tapering SIM images
PSFd = real(fftshift( ifft2(fftshift(OTFo.^3)) ));
PSFd = PSFd/max(max(PSFd));
PSFd = PSFd/sum(sum(PSFd));
h = 30;
PSFe = PSFd(wo-h+1:wo+h,wo-h+1:wo+h);

% edge tapering raw SIM images
S1aTnoisy = edgetaper(S1aTnoisy,PSFe);
S2aTnoisy = edgetaper(S2aTnoisy,PSFe);
S3aTnoisy = edgetaper(S3aTnoisy,PSFe);

fS1aTnoisy = fftshift(fft2(S1aTnoisy));
fS2aTnoisy = fftshift(fft2(S2aTnoisy));
fS3aTnoisy = fftshift(fft2(S3aTnoisy));

fDuplex2 = fS2aTnoisy - fS1aTnoisy;
fDuplex3 = fS3aTnoisy - fS1aTnoisy;

%% Optimizing the relative phases
Kai2Opt0 = @(phase0)Kai2Opt(phase0,fDuplex2,fDuplex3,OTFo);
options = optimset('LargeScale','off','Algorithm',...
	'active-set','MaxFunEvals',500,'MaxIter',500,'Display','notify');
phase0 = [pi/2 -pi/2]; % initial guess
[phaseShift,fval] = fminsearch(Kai2Opt0,phase0,options);
phaseShift*180/pi

%% Separating the three frequency components
phaseShift0 = 0;
[fDo,fDp,fDm] = SeparatedComponents2D(...
    phaseShift,phaseShift0,fS1aTnoisy,fS2aTnoisy,fS3aTnoisy);

%% for visual verification
%{
figure;
mesh(log(abs(fDo)))
figure;
mesh(log(abs(fDp)))
figure;
mesh(log(abs(fDm)))
% kk
%}