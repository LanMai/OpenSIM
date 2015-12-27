function Esum = OBJparaOpt(OBJpara0,fDIoTnoisy,OTFo)
% AIM: Determined Sum of Squared Errors (SSE) between `actual signal power' and
%   `approximated signal power'
% INPUT VARIABLES
%   OBJpara0: [Aobj Bobj], Object power parameters
%   fDIoTnoisy: FT of central frequency component
%   OTFo: system OTF
% OUTPUT VARIABLE
%   Esum: SSE between `actual signal spectrum' and `approximated signal spectrum'

w = size(OTFo,1);
wo = w/2;
x = linspace(0,w-1,w);
y = linspace(0,w-1,w);
[X,Y] = meshgrid(x,y);
Cv = (X-wo) + 1i*(Y-wo);
Ro = abs(Cv);
Ro(wo+1,wo+1) = 1; % to avoid nan

% approximated signal power calculation
Aobj = OBJpara0(1);
Bobj = OBJpara0(2);
OBJpower = Aobj*(Ro.^Bobj);
SIGpower = OBJpower.*OTFo;

% OTF cut-off frequency
Kotf = OTFedgeF(OTFo);

% range of frequency over which SSE is computed
Zloop = (Ro<0.75*Kotf).*(Ro>0.25*Kotf);

% frequency beyond which NoisePower estimate to be computed
NoiseFreq = Kotf + 20;

% NoisePower determination
Zo = Ro>NoiseFreq;
nNoise = fDIoTnoisy.*Zo;
NoisePower = sum(sum( nNoise.*conj(nNoise) ))./sum(sum(Zo));

% Noise free object power computation 
Fpower = fDIoTnoisy.*conj(fDIoTnoisy) - NoisePower;
fDIoTnoisy = sqrt(abs(Fpower));

% SSE computation
Error = fDIoTnoisy - SIGpower;
Esum = sum(sum((Error.^2./Ro).*Zloop));

