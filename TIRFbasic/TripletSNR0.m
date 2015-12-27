function [SIGao,SIGap2,SIGam2] = TripletSNR0(OBJpara,k2fa,OTFo,fDIp)

% AIM: To obtain signal spectrums corresponding to central and off-center
%       frequency components
% INPUT VARIABLES
%   OBJpara: object power parameters
%   k2fa: illumination frequency vector
%   OTFo: system OTF
%   fDIp: one of the off-center frequency component (utilized here only for
%           visual verification of computation)
% OUTPUT VARIABLES
%   SIGao: signal spectrum corresponding to central frequency component
%   SIGap2,SIGam2: signal spectrums corresponding to off-center frequency components


w = size(OTFo,1);
wo = w/2;
x = linspace(0,w-1,w);
y = linspace(0,w-1,w);
[X,Y] = meshgrid(x,y);
Cv = (X-wo) + 1i*(Y-wo);
Ro = abs(Cv);

% object spectrum parameters
Aobj = OBJpara(1);
Bobj = OBJpara(2);

% object spectrums (default)
kv = k2fa(2) + 1i*k2fa(1); % vector along illumination direction
Rp = abs(Cv-kv);
Rm = abs(Cv+kv);
OBJo = Aobj*(Ro.^Bobj);
OBJp = Aobj*(Rp.^Bobj);
OBJm = Aobj*(Rm.^Bobj);

OBJo(wo+1,wo+1) = 0.25*OBJo(wo+2,wo+1) + 0.25*OBJo(wo+1,wo+2)...
    + 0.25*OBJo(wo+0,wo+1) + 0.25*OBJo(wo+1,wo+0);

k3 = round(k2fa);
OBJp(wo+1+k3(1),wo+1+k3(2)) = 0.25*OBJp(wo+2+k3(1),wo+1+k3(2))...
	+ 0.25*OBJp(wo+1+k3(1),wo+2+k3(2))...
	+ 0.25*OBJp(wo+0+k3(1),wo+1+k3(2))...
	+ 0.25*OBJp(wo+1+k3(1),wo+0+k3(2));
OBJm(wo+1-k3(1),wo+1-k3(2)) = 0.25*OBJm(wo+2-k3(1),wo+1-k3(2))...
	+ 0.25*OBJm(wo+1-k3(1),wo+2-k3(2))...
	+ 0.25*OBJm(wo+0-k3(1),wo+1-k3(2))...
	+ 0.25*OBJm(wo+1-k3(1),wo+0-k3(2));

% signal spectrums
SIGao = OBJo.*OTFo;
SIGap = OBJp.*OTFo;
SIGam = OBJm.*OTFo;

SIGap2 = circshift(SIGap,-k3);
SIGam2 = circshift(SIGam,k3);

%% for visual verification
%{
pp = 3;
figure;
hold on
plot([0:w-1]-wo,log(abs(SIGap2(wo+pp,:))),'k--','LineWidth',3,'MarkerSize',6)
plot([0:w-1]-wo,log(abs(fDIp(wo+pp,:))),'mo-','LineWidth',2,'MarkerSize',6)
grid on
box on
%}