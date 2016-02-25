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
k2 = 75.23; % illumination freq
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

%% obtaining the noisy estimates of three frequency components
[fAo,fAp,fAm,kA]...
    = PCMseparateF(S1aTnoisy,S2aTnoisy,S3aTnoisy,OTFo);
% break
[fBo,fBp,fBm,kB]...
    = PCMseparateF(S1bTnoisy,S2bTnoisy,S3bTnoisy,OTFo);
[fCo,fCp,fCm,kC]...
    = PCMseparateF(S1cTnoisy,S2cTnoisy,S3cTnoisy,OTFo);

% averaging the central frequency components
fCent = (fAo + fBo + fCo)/3;

% Object power parameters determination
OBJparaA = OBJpowerPara(fCent,OTFo);


%% Wiener Filtering the noisy frequency components
[fAof,fApf,fAmf,Nao,Nap,Nam,Ma,DoubleMatSize]...
    = PCMfilteringF(fAo,fAp,fAm,OTFo,OBJparaA,kA);
[fBof,fBpf,fBmf,Nbo,Nbp,Nbm,Mb,~]...
    = PCMfilteringF(fBo,fBp,fBm,OTFo,OBJparaA,kB);
[fCof,fCpf,fCmf,Nco,Ncp,Ncm,Mc,~]...
    = PCMfilteringF(fCo,fCp,fCm,OTFo,OBJparaA,kC);

%% doubling Fourier domain size if necessary
OTFo = OTFdoubling(OTFo,DoubleMatSize);

%% merging all 9 frequency components using generalized Wiener Filter
[Fsum,Fperi,Fcent] = MergingHeptaletsF(fAof,fApf,fAmf,...
    fBof,fBpf,fBmf,fCof,fCpf,fCmf,...
    Ma,Mb,Mc,Nao,Nap,Nam,Nbo,Nbp,Nbm,...
    Nco,Ncp,Ncm,kA,kB,kC,OBJparaA,OTFo);

% Plotting SIM results
SIMplot(Fsum,Fperi,Fcent,OTFo,kA,kB,kC,S1aTnoisy);
%break

%% recontructed SIM images
Dcent = real( ifft2(fftshift(Fcent)) );
Dsum = real( ifft2(fftshift(Fsum)) );
Dperi = real( ifft2(fftshift(Fperi)) );

%{
figure;
mesh(log(abs(Fsum)))
figure;
mesh(log(abs(Fperi)))
%}

w = size(OTFo,1);
h = 30; % pixel width to be trimed of the image edge to trim of artifacts
figure;
imshow(Dcent(h+1:t-h,h+1:t-h),[])
title('Wiener Filtered wide-field')
figure;
imshow(Dsum(h+1:t-h,h+1:t-h),[])
title('SIM image')
figure;
imshow(Dperi(h+1:t-h,h+1:t-h),[])
title('SIM image (using only off-center frequency components)')

%% appodizing the merged frequency components
Index = 0.4;
Kotf = OTFedgeF(OTFo);
[FsumA] = ApodizationFunction(Fsum,kA,kB,kC,Kotf,Index);
[FperiA] = ApodizationFunction(Fperi,kA,kB,kC,Kotf,Index);
DsumA = real( ifft2(fftshift(FsumA)) );
DperiA = real( ifft2(fftshift(FperiA)) );

figure;
imshow(DsumA(h+1:t-h,h+1:t-h),[])
title('appodized SIM image')
figure;
imshow(DperiA(h+1:t-h,h+1:t-h),[])
title('appodized SIM image (using only off-center frequency components)')
%break
%% for writing the image files into *.tiff files
%{
h = 1*30;
imwrite(im2uint16(mat2gray(Dcent(h+1:t-h,h+1:t-h))),['Dcent.tif'],'tiff');
imwrite(im2uint16(mat2gray(Dsum(h+1:t-h,h+1:t-h))),['Dsum.tif'],'tiff');
imwrite(im2uint16(mat2gray(Dperi(h+1:t-h,h+1:t-h))),['Dperi.tif'],'tiff');
imwrite(im2uint16(mat2gray(DsumA(h+1:t-h,h+1:t-h))),['DsumA.tif'],'tiff');
imwrite(im2uint16(mat2gray(DperiA(h+1:t-h,h+1:t-h))),['DperiA.tif'],'tiff');
%}

%% ploting FT of images
fDIo = fftshift(fft2(DIo));
fDIoT = fftshift(fft2(DIoT));
fDIoTnoisy = fftshift(fft2(DIoTnoisy));
fS1aTnoisy = fftshift(fft2(S1aTnoisy));

p = 10;
minL1 = min(min( abs(fDIo).^(1/p) ));
minL2 = min(min( abs(fDIoT).^(1/p) ));
minL3 = min(min( abs(fDIoTnoisy).^(1/p) ));
minL4 = min(min( abs(fS1aTnoisy).^(1/p) ));
minL5 = min(min( abs(Fcent).^(1/p) ));
minL6 = min(min( abs(Fsum).^(1/p) ));
minL7 = min(min( abs(Fperi).^(1/p) ));
maxL1 = max(max( abs(fDIo).^(1/p) ));
maxL2 = max(max( abs(fDIoT).^(1/p) ));
maxL3 = max(max( abs(fDIoTnoisy).^(1/p) ));
maxL4 = max(max( abs(fS1aTnoisy).^(1/p) ));
maxL5 = max(max( abs(Fcent).^(1/p) ));
maxL6 = max(max( abs(Fsum).^(1/p) ));
maxL7 = max(max( abs(Fperi).^(1/p) ));
minL = min([minL1,minL2,minL3,minL4,minL5,minL6,minL7]);
maxL = max([maxL1,maxL2,maxL3,maxL4,maxL5,maxL6,maxL7]);

figure;
imshow(abs(fDIo).^(1/p),[minL maxL])
title('fDIo')
figure;
imshow(abs(fDIoT).^(1/p),[minL maxL])
title('fDIoT')
figure;
imshow(abs(fDIoTnoisy).^(1/p),[minL maxL])
title('fDIoTnoisy')
figure;
imshow(abs(fS1aTnoisy).^(1/p),[minL maxL])
title('fS1aTnoisy')

figure;
imshow(abs(Fcent).^(1/p),[minL maxL])
title('Fcent')
figure;
imshow(abs(Fsum).^(1/p),[minL maxL])
title('Fsum')
figure;
imshow(abs(Fperi).^(1/p),[minL maxL])
title('Fperi')
figure;
imshow(abs(FsumA).^(1/p),[minL maxL])
title('FsumA')
figure;
imshow(abs(FperiA).^(1/p),[minL maxL])
title('FperiA')

%% for writing the image files into *.tiff files
%{
imwrite(im2uint16(mat2gray(abs(Fcent).^(1/p))),['Fcent.tif'],'tiff');
imwrite(im2uint16(mat2gray(abs(Fsum).^(1/p))),['Fsum.tif'],'tiff');
imwrite(im2uint16(mat2gray(abs(Fperi).^(1/p))),['Fperi.tif'],'tiff');
imwrite(im2uint16(mat2gray(abs(FsumA).^(1/p))),['FsumA.tif'],'tiff');
imwrite(im2uint16(mat2gray(abs(FperiA).^(1/p))),['FperiA.tif'],'tiff');
%}
