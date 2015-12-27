function OTFo = OTFpost(OTFo)
%% AIM: To flaten out the peripheral rippling frequencies

OTFo = OTFo./max(OTFo(:));

% OTF cut-off frequency
Kotf = OTFedgeF(OTFo);

w = size(OTFo,1);
wo = w/2;
x = linspace(0,w-1,w);
y = linspace(0,w-1,w);
[X,Y] = meshgrid(x,y);
Cv = (X-wo) + 1i*(Y-wo);
Ro = abs(Cv);

Zmask = (Ro>(Kotf+20)).*(Ro<wo);
OTFperi = sum(sum(OTFo.*Zmask))./sum(sum(Zmask));
Zperi = (Ro>(wo-10));
OTFo = OTFo.*(1-Zperi) + OTFperi.*Zperi;