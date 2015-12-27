function [yy0,OTF2dc] = PsfOtf(w,scale)

% AIM: To generate PSF and OTF using Bessel function
% INPUT VARIABLES
%   w: image size
%   scale: a parameter used to adjust PSF/OTF width
% OUTPUT VRAIBLES
%   yyo: system PSF
%   OTF2dc: system OTF

x = linspace(0,w-1,w);
y = linspace(0,w-1,w);
[X,Y] = meshgrid(x,y);

%% Generation of the PSF with Besselj.
% scale = 0.3;
R=sqrt(min(X,abs(X-w)).^2+min(Y,abs(Y-w)).^2);
yy=abs(2*besselj(1,scale*R+eps,1)./(scale*R+eps)).^2; %0.5 is introduced to make PSF wider
yy0 = fftshift(yy);
% figure;
% mesh(X,Y,yy0)
% figure;
% imshow(yy0,[])
% axis square

%Generate 2D OTF.
OTF2d=fft2(yy);
OTF2dmax = max(max(abs(OTF2d)));
OTF2d = OTF2d./OTF2dmax;
OTF2dc = abs(fftshift(OTF2d));
% figure;
% mesh(abs(OTF2dc))
% title('OTF of the system')
% axis square