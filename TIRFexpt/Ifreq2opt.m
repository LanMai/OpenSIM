function [CCop] = Ifreq2opt(k2fa,fAo,fAp,OTFo)
% Aim: compute cross-correlation between central and off-central frequency
% components 
% INPUT VARIABLES:
%   k2fa: nearest pixel approximation to illumination frequency vector
%   fAo: central frequency component
%   fAp: off-center frequency component
%   OTFo: system OTF

w = size(OTFo,1);
wo = w/2;
x = linspace(0,w-1,w);
y = linspace(0,w-1,w);
[X,Y] = meshgrid(x,y);
Cv = (X-wo) + 1i*(Y-wo);
Ro = abs(Cv);

% suppressing residual low frequency peaks of off-center components using
% G is a notch-filter (determined heuristically)
G = 1 - exp(-0.05*Ro.^1.2);
fAo = fAo.*G;
fAp = fAp.*G;

fAoT = fAo.*conj(OTFo);
fApT = fAp.*conj(OTFo);

Kotf = OTFedgeF(OTFo);
DoubleMatSize = 0; 
if ( 2*Kotf > wo )
	DoubleMatSize = 1; % 1 for doubling fourier domain size, 0 for keeping it unchanged
end

if ( DoubleMatSize>0 )
    t = 2*w;
    fAoT_temp = zeros(t,t);
    fApT_temp = zeros(t,t);
    fAoT_temp(wo+1:w+wo,wo+1:w+wo) = fAoT;
    fApT_temp(wo+1:w+wo,wo+1:w+wo) = fApT;
    clear fAoT fApT
    fAoT = fAoT_temp;
    fApT = fApT_temp;
    clear fAoT_temp fApT_temp
else
    t = w;
end
to = t/2;
u = linspace(0,t-1,t);
v = linspace(0,t-1,t);
[U,V] = meshgrid(u,v);
ApT = exp( +1i.*2*pi*( k2fa(2)/t.*(U-to)+k2fa(1)/t.*(V-to) ) ).*ifft2(fApT);
fApT0 = fft2( ApT );

mA = sum(sum( fAoT.*conj(fApT0) ));
mA = mA./sum(sum( fApT0.*conj(fApT0) ));

CCop = -abs(mA);