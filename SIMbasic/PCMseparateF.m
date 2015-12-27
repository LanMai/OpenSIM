function [fDo,fDp,fDm,kA]...
    = PCMseparateF(S1aTnoisy,S2aTnoisy,S3aTnoisy,OTFo)
% AIM: obtaining the noisy estimates of three frequency components
% INPUT VARIABLES
%   S1aTnoisy,S2aTnoisy,S3aTnoisy: 3 raw SIM images with identical 
%                                   illumination pattern orientation 
%                                   but different phase shifts
%   OTFo: system OTF
% OUTPUT VARIABLES
%   fDo,fDp,fDm: noisy estimates of separated frequency components
%   kA: (averaged) illumination frequency vector


w = size(S1aTnoisy,1);
wo = w/2;

%% Determination of illumination frequency vectors
[k1a] = IlluminationFreqF(S1aTnoisy,OTFo);
[k2a] = IlluminationFreqF(S2aTnoisy,OTFo);
[k3a] = IlluminationFreqF(S3aTnoisy,OTFo);

% mean illumination frequency vector
kA = (k1a + k2a + k3a)/3;

%% determination of illumination phase shifts
[phase1A] = IlluminationPhaseF(S1aTnoisy,kA);
[phase2A] = IlluminationPhaseF(S2aTnoisy,kA);
[phase3A] = IlluminationPhaseF(S3aTnoisy,kA);

% for display in command window
phaseA = [phase1A; phase2A; phase3A];
phaseA*180/pi

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

%% Separating the three frequency components
phaseShift0 = phase1A;
phaseShift = [phase2A; phase3A];
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
%}
