function Mm = ModulationFactor(fDp,k1a,OBJparaA,OTFo)
% AIM: Determination of modulation factor
% INPUT VARIABLES
%   fDp: off-center frequency component
%   k1a: illumination frequency vector
%   OBJparaA: Object power parameters
%   OTFo: system OTF
% OUTPUT VARIABLE
%   Mm: modulation factor

w = size(OTFo,1);
wo = w/2;
x = linspace(0,w-1,w);
y = linspace(0,w-1,w);
[X,Y] = meshgrid(x,y);
Cv = (X-wo) + 1i*(Y-wo);
Ro = abs(Cv);

% magnitude of illumination vector
k2 = sqrt(k1a*k1a');

% vector along illumination direction
kv = k1a(2) + 1i*k1a(1); 
Rp = abs(Cv+kv);

% Object power parameters
Aobj = OBJparaA(1);
Bobj = OBJparaA(2);

% Object spectrum
OBJp = Aobj*(Rp+0).^Bobj;

% illumination vector rounded to nearest pixel
k3 = -round(k1a);

OBJp(wo+1+k3(1),wo+1+k3(2)) = 0.25*OBJp(wo+2+k3(1),wo+1+k3(2))...
	+ 0.25*OBJp(wo+1+k3(1),wo+2+k3(2))...
	+ 0.25*OBJp(wo+0+k3(1),wo+1+k3(2))...
	+ 0.25*OBJp(wo+1+k3(1),wo+0+k3(2));

% signal spectrum
SIGap = OBJp.*OTFo;

% figure;
% mesh(log(abs(OBJp)))
% figure;
% mesh(log(abs(fDp)))

% OTF cut-off frequency
Kotf = OTFedgeF(OTFo);

% frequency beyond which NoisePower estimate to be computed
NoiseFreq = Kotf + 20;

% NoisePower determination
Zo = Ro>NoiseFreq;
nNoise = fDp.*Zo;
NoisePower = sum(sum( nNoise.*conj(nNoise) ))./sum(sum(Zo));

% Noise free object power computation
Fpower = fDp.*conj(fDp) - NoisePower;
fDp = sqrt(abs(Fpower));

% frequency range over which signal power matching is done to estimate
% modulation factor
Zmask = (Ro > 0.2*k2).*(Ro < 0.8*k2).*(Rp > 0.2*k2);

% least square approximation for modulation factor
Mm = sum(sum(SIGap.*abs(fDp).*Zmask));
Mm = Mm./sum(sum(SIGap.^2.*Zmask))

%% for visual verification
%{
pp = 3;
figure;
hold on
plot(x-wo,log(abs(fDp(wo+pp,:))),'o-','LineWidth',2,'MarkerSize',6)
plot(x-wo,log(abs(SIGap(wo+pp,:))),'go--','LineWidth',2,'MarkerSize',4)
plot(x-wo,log(abs(Mm*SIGap(wo+pp,:))),'r*--','LineWidth',2,'MarkerSize',4)
grid on
box on
kk
%}



