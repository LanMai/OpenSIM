function [OTFap_mask] = OTFmaskShifted(k2fa,Kotf,w)
% OTFap_mask = complement of shifted OTF

wo = w/2;
x = linspace(0,w-1,w);
y = linspace(0,w-1,w);
[X,Y] = meshgrid(x,y);

Iy = k2fa(1,1) + (wo+1);
Ix = k2fa(1,2) + (wo+1);

Ro = sqrt( (X+1-Ix).^2 + (Y+1-Iy).^2 );
OTFap_mask = Ro > Kotf;