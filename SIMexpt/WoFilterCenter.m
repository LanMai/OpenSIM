function [FiSMaof,NoisePower] = WoFilterCenter(FiSMao,OTFo,co,OBJpara,SFo)
% Aim: Wiener Filtering the central frequency component
% INPUT VARIABLES
%   FiSMao: noisy central frequency component
%   OTFo: system OTF
%   co: Wiener filter constant [=1, for minimum RMS estimate]
%   SFo: scaling factor (not significant here, so set to 1)
% OUTPUT VARIABLES
%   FiSMaof: Wiener Filtered estimate of FiSMao
%   NoisePower: avg. noise power in FiSMao

w = size(FiSMao,1);
wo = w/2;
x = linspace(0,w-1,w);
y = linspace(0,w-1,w);
[X,Y] = meshgrid(x,y);
Ro = sqrt( (X-wo).^2 + (Y-wo).^2 );

OTFpower = OTFo.*conj(OTFo);

% OTF cut-off frequency
Kotf = OTFedgeF(OTFo);

% frequency beyond which NoisePower estimate to be computed
NoiseFreq = Kotf + 20;

% NoisePower determination
Zo = Ro>NoiseFreq;
nNoise = FiSMao.*Zo;
NoisePower = sum(sum( nNoise.*conj(nNoise) ))./sum(sum(Zo));

% Object Power determination
Aobj = OBJpara(1);
Bobj = OBJpara(2);
Ro(wo+1,wo+1) = 1;
OBJpower = Aobj*(Ro.^Bobj);
OBJpower = OBJpower.^2;

%% Wiener Filtering
FiSMaof = FiSMao.*(SFo.*conj(OTFo)./NoisePower)./((SFo.^2).*OTFpower./NoisePower + co./OBJpower);

%% for cross-checking filtered estimate visually
%{
WFilter = (SFo.^2.*OTFpower./NoisePower)./((SFo.^2).*OTFpower./NoisePower + co./OBJpower);
FiSMao1 = FiSMaof.*OTFo;
pp = 3;
figure;
hold on
plot(x-wo,log(abs(FiSMaof(wo+pp,:))),'go--','LineWidth',2,'MarkerSize',6)
plot(x-wo,log(abs(FiSMao(wo+pp,:))),'o--','LineWidth',2,'MarkerSize',6)
plot(x-wo,log(abs(FiSMao1(wo+pp,:))),'r*--','LineWidth',2,'MarkerSize',4)
plot(x-wo,10*WFilter(wo+pp,:),'k*--','LineWidth',2,'MarkerSize',4)
title('WoFilter')
legend('FiSMaof','FiSMao','FiSMao1')
grid on
box on
%}