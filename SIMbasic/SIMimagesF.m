function [S1aTnoisy S2aTnoisy S3aTnoisy ...
          S1bTnoisy S2bTnoisy S3bTnoisy ...
          S1cTnoisy S2cTnoisy S3cTnoisy ...
          DIoTnoisy DIoT] = SIMimagesF(k2,...
          DIo,PSFo,OTFo,ModFac,NoiseLevel,UsePSF)

% AIM: to generate raw sim images
% INPUT VARIABLES
%   k2: illumination frequency
%   DIo: specimen image
%   PSFo: system PSF
%   OTFo: system OTF
%   UsePSF: 1 (to blur SIM images by convloving with PSF)
%           0 (to blur SIM images by truncating its fourier content beyond OTF)
%   NoiseLevel: percentage noise level for generating gaussian noise
% OUTPUT VARIABLES
%   [S1aTnoisy S2aTnoisy S3aTnoisy
%    S1bTnoisy S2bTnoisy S3bTnoisy
%    S1cTnoisy S2cTnoisy S3cTnoisy]: nine raw sim images
%   DIoTnoisy: noisy wide field image
%   DIoT: noise-free wide field image
      
w = size(DIo,1);
wo = w/2;
x = linspace(0,w-1,w);
y = linspace(0,w-1,w);
[X,Y] = meshgrid(x,y);

%% illunination phase shifts along the three directions
p0Ao = 0*pi/3;
p0Ap = 2*pi/3;
p0Am = 4*pi/3;
p0Bo = 0*pi/3;
p0Bp = 2*pi/3;
p0Bm = 4*pi/3;
p0Co = 0*pi/3;
p0Cp = 2*pi/3;
p0Cm = 4*pi/3;

%% Illuminating patterns
alpha = 0*pi/6; 
% orientation direction of illumination patterns
thetaA = 0*pi/3 + alpha; 
thetaB = 1*pi/3 + alpha;
thetaC = 2*pi/3 + alpha;
% illumination frequency vectors
k2a = (k2/w).*[cos(thetaA) sin(thetaA)];
k2b = (k2/w).*[cos(thetaB) sin(thetaB)];
k2c = (k2/w).*[cos(thetaC) sin(thetaC)];
% -------------------------------------------------------
% mean illumination intensity
mA = 0.5; 
mB = 0.5;
mC = 0.5;
% amplitude of illumination intensity above mean
aA = 0.5*ModFac; 
aB = 0.5*ModFac;
aC = 0.5*ModFac;

% random phase shift errors
NN = 1*(0.5-rand(9,1))*pi/18;

% illunination phase shifts with random errors
psAo = p0Ao + NN(1,1);
psAp = p0Ap + NN(2,1);
psAm = p0Am + NN(3,1);
psBo = p0Bo + NN(4,1);
psBp = p0Bp + NN(5,1);
psBm = p0Bm + NN(6,1);
psCo = p0Co + NN(7,1);
psCp = p0Cp + NN(8,1);
psCm = p0Cm + NN(9,1);

% illunination patterns
sAo = mA + aA*cos(2*pi*(k2a(1,1).*(X-wo)+k2a(1,2).*(Y-wo))+psAo); % illuminated signal (0 phase)
sAp = mA + aA*cos(2*pi*(k2a(1,1).*(X-wo)+k2a(1,2).*(Y-wo))+psAp); % illuminated signal (+ phase)
sAm = mA + aA*cos(2*pi*(k2a(1,1).*(X-wo)+k2a(1,2).*(Y-wo))+psAm); % illuminated signal (- phase)
sBo = mB + aB*cos(2*pi*(k2b(1,1).*(X-wo)+k2b(1,2).*(Y-wo))+psBo); % illuminated signal (0 phase)
sBp = mB + aB*cos(2*pi*(k2b(1,1).*(X-wo)+k2b(1,2).*(Y-wo))+psBp); % illuminated signal (+ phase)
sBm = mB + aB*cos(2*pi*(k2b(1,1).*(X-wo)+k2b(1,2).*(Y-wo))+psBm); % illuminated signal (- phase)
sCo = mC + aC*cos(2*pi*(k2c(1,1).*(X-wo)+k2c(1,2).*(Y-wo))+psCo); % illuminated signal (0 phase)
sCp = mC + aC*cos(2*pi*(k2c(1,1).*(X-wo)+k2c(1,2).*(Y-wo))+psCp); % illuminated signal (+ phase)
sCm = mC + aC*cos(2*pi*(k2c(1,1).*(X-wo)+k2c(1,2).*(Y-wo))+psCm); % illuminated signal (- phase)

%{
fftemp = fftshift(fft2(sAo));
figure;
hold on
plot(x-wo,abs(OTFo(wo+1,:)),'--','LineWidth',2,'MarkerSize',6)
plot(x-wo,0.5*abs(fftemp(wo+1,:))./max(max(abs(fftemp))),'r--','LineWidth',2,'MarkerSize',6)
legend('OTFo','illumination spectrum')
grid on
box on
%}

% figure;
% imshow(sAo,[ ])
% kk

%% superposed Objects
s1a = DIo.*sAo; % superposed signal (0 phase)
s2a = DIo.*sAp; % superposed signal (+ phase)
s3a = DIo.*sAm; % superposed signal (- phase)
s1b = DIo.*sBo; 
s2b = DIo.*sBp; 
s3b = DIo.*sBm; 
s1c = DIo.*sCo; 
s2c = DIo.*sCp; 
s3c = DIo.*sCm;
%{
figure;
subplot(3,3,1)
imshow(s1a,[ ])
title('s1a')
subplot(3,3,4)
imshow(s2a,[ ])
title('s2a')
subplot(3,3,7)
imshow(s3a,[ ])
title('s3a')
subplot(3,3,2)
imshow(s1b,[ ])
title('s1b')
subplot(3,3,5)
imshow(s2b,[ ])
title('s2b')
subplot(3,3,8)
imshow(s3b,[ ])
title('s3b')
subplot(3,3,3)
imshow(s1c,[ ])
title('s1c')
subplot(3,3,6)
imshow(s2c,[ ])
title('s2c')
subplot(3,3,9)
imshow(s3c,[ ])
title('s3c')
%}


%% superposed (noise-free) Images
PSFsum = sum(sum(PSFo));
if ( UsePSF == 1 )
    DIoT = conv2(DIo,PSFo,'same')./PSFsum;
    S1aT = conv2(s1a,PSFo,'same')./PSFsum;
    S2aT = conv2(s2a,PSFo,'same')./PSFsum;
    S3aT = conv2(s3a,PSFo,'same')./PSFsum;
    S1bT = conv2(s1b,PSFo,'same')./PSFsum;
    S2bT = conv2(s2b,PSFo,'same')./PSFsum;
    S3bT = conv2(s3b,PSFo,'same')./PSFsum;
    S1cT = conv2(s1c,PSFo,'same')./PSFsum;
    S2cT = conv2(s2c,PSFo,'same')./PSFsum;
    S3cT = conv2(s3c,PSFo,'same')./PSFsum;
else
    DIoT = ifft2( fft2(DIo).*fftshift(OTFo) );
    S1aT = ifft2( fft2(s1a).*fftshift(OTFo) );
    S2aT = ifft2( fft2(s2a).*fftshift(OTFo) );
    S3aT = ifft2( fft2(s3a).*fftshift(OTFo) );
    S1bT = ifft2( fft2(s1b).*fftshift(OTFo) );
    S2bT = ifft2( fft2(s2b).*fftshift(OTFo) );
    S3bT = ifft2( fft2(s3b).*fftshift(OTFo) );
    S1cT = ifft2( fft2(s1c).*fftshift(OTFo) );
    S2cT = ifft2( fft2(s2c).*fftshift(OTFo) );
    S3cT = ifft2( fft2(s3c).*fftshift(OTFo) );
    
    DIoT = real(DIoT);
    S1aT = real(S1aT);
    S2aT = real(S2aT);
    S3aT = real(S3aT);
    S1bT = real(S1bT);
    S2bT = real(S2bT);
    S3bT = real(S3bT);
    S1cT = real(S1cT);
    S2cT = real(S2cT);
    S3cT = real(S3cT);    
end
%{
figure;
subplot(3,3,1)
imshow(S1aT,[ ])
title('S1aT')
subplot(3,3,4)
imshow(S2aT,[ ])
title('S2aT')
subplot(3,3,7)
imshow(S3aT,[ ])
title('S3aT')
subplot(3,3,2)
imshow(S1bT,[ ])
title('S1bT')
subplot(3,3,5)
imshow(S2bT,[ ])
title('S2bT')
subplot(3,3,8)
imshow(S3bT,[ ])
title('S3bT')
subplot(3,3,3)
imshow(S1cT,[ ])
title('S1cT')
subplot(3,3,6)
imshow(S2cT,[ ])
title('S2cT')
subplot(3,3,9)
imshow(S3cT,[ ])
title('S3cT')
%}


%% Gaussian noise generation
aNoise = NoiseLevel/100; % corresponds to 10% noise
% SNR = 1/aNoise
% SNRdb = 20*log10(1/aNoise)
nDIoT = random('norm', 0, aNoise*std2(DIoT), w , w);
nS1aT = random('norm', 0, aNoise*std2(S1aT), w , w);
nS2aT = random('norm', 0, aNoise*std2(S2aT), w , w);
nS3aT = random('norm', 0, aNoise*std2(S3aT), w , w);
nS1bT = random('norm', 0, aNoise*std2(S2bT), w , w);
nS2bT = random('norm', 0, aNoise*std2(S2bT), w , w);
nS3bT = random('norm', 0, aNoise*std2(S2bT), w , w);
nS1cT = random('norm', 0, aNoise*std2(S3cT), w , w);
nS2cT = random('norm', 0, aNoise*std2(S3cT), w , w);
nS3cT = random('norm', 0, aNoise*std2(S3cT), w , w);

%% noise added raw SIM images
NoiseFrac = 1; %may be set to 0 to avoid noise addition
DIoTnoisy = DIoT + NoiseFrac*nDIoT;
S1aTnoisy = S1aT + NoiseFrac*nS1aT;
S2aTnoisy = S2aT + NoiseFrac*nS2aT;
S3aTnoisy = S3aT + NoiseFrac*nS3aT;
S1bTnoisy = S1bT + NoiseFrac*nS1bT;
S2bTnoisy = S2bT + NoiseFrac*nS2bT;
S3bTnoisy = S3bT + NoiseFrac*nS3bT;
S1cTnoisy = S1cT + NoiseFrac*nS1cT;
S2cTnoisy = S2cT + NoiseFrac*nS2cT;
S3cTnoisy = S3cT + NoiseFrac*nS3cT;
%{
figure;
subplot(3,3,1)
imshow(S1aTnoisy,[ ])
title('S1aTnoisy')
subplot(3,3,4)
imshow(S2aTnoisy,[ ])
title('S2aTnoisy')
subplot(3,3,7)
imshow(S3aTnoisy,[ ])
title('S3aTnoisy')
subplot(3,3,2)
imshow(S1bTnoisy,[ ])
title('S1bTnoisy')
subplot(3,3,5)
imshow(S2bTnoisy,[ ])
title('S2bTnoisy')
subplot(3,3,8)
imshow(S3bTnoisy,[ ])
title('S3bTnoisy')
subplot(3,3,3)
imshow(S1cTnoisy,[ ])
title('S1cTnoisy')
subplot(3,3,6)
imshow(S2cTnoisy,[ ])
title('S2cTnoisy')
subplot(3,3,9)
imshow(S3cTnoisy,[ ])
title('S3cTnoisy')
%}

%{
figure;
imshow(S1aT,[ ])
title('S1aT')
figure;
imshow(S1aTnoisy,[ ])
title('S1aTnoisy')

figure;
hold on
plot(x-wo,S1aT(wo+1,:),'o-','LineWidth',2,'MarkerSize',6)
plot(x-wo,S1aTnoisy(wo+1,:),'r*--','LineWidth',2,'MarkerSize',6)
grid on
box on
%}