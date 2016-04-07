% Copyright (c) <2016-2018> <Amit Lal, Peng Xi>
% 
% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the "Software"), 
% to deal in the Software without restriction, including without limitation 
% the rights to use, copy, modify, merge, publish, distribute, sublicense, 
% and/or sell copies of the Software, and to permit persons to whom the 
% Software is furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
% DEALINGS IN THE SOFTWARE.

clear all
close all
clc

w = 512;
wo = w/2;
x = linspace(0,w-1,w);
y = linspace(0,w-1,w);
[X,Y] = meshgrid(x,y);

%% Generation of the PSF with Besselj.
scale = 0.63; % used to adjust PSF/OTF width
[PSFo,OTFo] = PsfOtf(w,scale); 
%{
figure;
mesh(PSFo)
title('system PSF')

figure;
mesh(OTFo)
title('system OTF')
%}

%% Reading input file
Io1 = imread('testpat.tiff');
Io = Io1(257:768,257:768); % selecting central 512x512 of image
DIo = double(Io);
figure;
imshow(DIo,[])


%% Generating raw SIM Images
k2 = 120.23; % illumination freq
ModFac = 0.8; % modulation factor
NoiseLevel = 10; % in percentage
UsePSF = 0; % 1(to blur using PSF), 0(to blur using OTF)
[S1aTnoisy S2aTnoisy S3aTnoisy ...
 S1bTnoisy S2bTnoisy S3bTnoisy ...
 S1cTnoisy S2cTnoisy S3cTnoisy ...
 DIoTnoisy DIoT] = SIMimagesF(k2,DIo,PSFo,OTFo,ModFac,NoiseLevel,UsePSF);
%{
figure;
imshow(S1aTnoisy,[])
title('Raw SIM image')
figure;
imshow(S1bTnoisy,[])
title('Raw SIM image')
figure;
imshow(S1cTnoisy,[])
title('Raw SIM image')
figure;
imshow(DIoTnoisy,[])
title('Noisy Wide-Field image')
figure;
imshow(DIoT,[])
title('Noise-free Wide-Field image')
break
%}

% tentative values to define search region for illumination frequency
% determination
k2o = 125; % tentative illumination frequency 
thetaA = 0*pi/3; % orientations of structured illumination
thetaB = 1*pi/3;
thetaC = 2*pi/3;

%% obtaining the noisy estimates of three frequency components
[fAo,fAp,fAm]...
    = PCMseparateF(S1aTnoisy,S2aTnoisy,S3aTnoisy,OTFo);
kA = IlluminationFreqTIRF(fAo,fAp,OTFo,k2o,thetaA);
[fBo,fBp,fBm]...
    = PCMseparateF(S1bTnoisy,S2bTnoisy,S3bTnoisy,OTFo);
kB = IlluminationFreqTIRF(fBo,fBp,OTFo,k2o,thetaB);
[fCo,fCp,fCm]...
    = PCMseparateF(S1cTnoisy,S2cTnoisy,S3cTnoisy,OTFo);
kC = IlluminationFreqTIRF(fCo,fCp,OTFo,k2o,thetaC);

% averaging the central frequency components
fCent = (fAo + fBo + fCo)/3;

% Object power parameters determination
OBJparaA = OBJpowerPara(fCent,OTFo);

% break
%% Wiener Filtering the noisy frequency components
[fAof,fApf,fAmf,Nao,Nap,Nam,Ma,DoubleMatSize]...
    = PCMfilteringF(fAo,fAp,fAm,OTFo,OBJparaA,kA);
[fBof,fBpf,fBmf,Nbo,Nbp,Nbm,Mb,~]...
    = PCMfilteringF(fAo,fBp,fBm,OTFo,OBJparaA,kB);
[fCof,fCpf,fCmf,Nco,Ncp,Ncm,Mc,~]...
    = PCMfilteringF(fAo,fCp,fCm,OTFo,OBJparaA,kC);

%% doubling Fourier domain size if necessary
OTFo = OTFdoubling(OTFo,DoubleMatSize);

%% merging all 9 frequency components using generalized Wiener Filter
[Fsum,Fperi,Fcent] = MergingHeptaletsF(fAof,fApf,fAmf,...
    fBof,fBpf,fBmf,fCof,fCpf,fCmf,...
    Ma,Mb,Mc,Nao,Nap,Nam,Nbo,Nbp,Nbm,...
    Nco,Ncp,Ncm,kA,kB,kC,OBJparaA,OTFo);

% Plotting SIM results
SIMplot(Fsum,Fperi,Fcent,OTFo,kA,kB,kC,S1aTnoisy);

