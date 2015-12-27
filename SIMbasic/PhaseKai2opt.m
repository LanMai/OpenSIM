function [CCop] = PhaseKai2opt(k2fa,fS1aTnoisy,OTFo,OPT)
% Aim: Compute autocorrelation of FT of raw SIM images
%   k2fa: illumination frequency vector
%   fS1aTnoisy: FT of raw SIM image
%   OTFo: system OTF
%   OPT: acronym for `OPTIMIZE'; to be set to 1 when this function is used
%       for optimization, or else to 0
%   CCop: autocorrelation of fS1aTnoisy

w = size(fS1aTnoisy,1);
wo = w/2;

fS1aTnoisy = fS1aTnoisy.*(1-1*OTFo.^10);
fS1aT = fS1aTnoisy.*conj(OTFo);

Kotf = OTFedgeF(OTFo);
DoubleMatSize = 0; 
if ( 2*Kotf > wo )
	DoubleMatSize = 1; % 1 for doubling fourier domain size, 0 for keeping it unchanged
end

if ( DoubleMatSize>0 )
    t = 2*w;
    fS1aT_temp = zeros(t,t);
    fS1aT_temp(wo+1:w+wo,wo+1:w+wo) = fS1aT;
    clear fS1aT
    fS1aT = fS1aT_temp;
    clear fS1aT_temp
else
    t = w;
end
to = t/2;
u = linspace(0,t-1,t);
v = linspace(0,t-1,t);
[U,V] = meshgrid(u,v);
S1aT = exp( -1i.*2*pi*( k2fa(2)/t.*(U-to)+k2fa(1)/t.*(V-to) ) ).*ifft2(fS1aT);
fS1aT0 = fft2( S1aT );

mA = sum(sum( fS1aT.*conj(fS1aT0) ));
mA = mA./sum(sum( fS1aT0.*conj(fS1aT0) ));

if (OPT > 0)
    CCop = -abs(mA);
else
    CCop = mA;    
end