function [maxK2,Ix,Iy] = ApproxFreqDuplex(FiSMap,Kotf)
% AIM: approx. illumination frequency vector determination
% INPUT VARIABLES
%   FiSMap: FT of raw SIM image
%   Kotf: OTF cut-off frequency
% OUTPUT VARIABLES
%   maxK2: illumination frequency vector (approx)
%   Ix,Iy: coordinates of illumination frequency peaks

FiSMap = abs(FiSMap);

w = size(FiSMap,1);
wo = w/2;
x = linspace(0,w-1,w);
y = linspace(0,w-1,w);
[X,Y] = meshgrid(x,y);

Ro = sqrt( (X-wo).^2 + (Y-wo).^2 );
Z0 = Ro > round(0.5*Kotf);
Z1 = X > wo;

FiSMap = FiSMap.*Z0.*Z1;
dumY = max( FiSMap,[],1 );
[dummy Iy] = max(dumY);
dumX = max( FiSMap,[],2 );
[dummy Ix] = max(dumX);

maxK2 = [Ix-(wo+1) Iy-(wo+1)];
