function [FiSMaof,NoisePower] = WoFilterSideLobe(FiSMao,OTFo,co,OBJsideP,SFo)
% Aim: Wiener Filtering the off-center frequency component
% INPUT VARIABLES
%   FiSMao: noisy off-center frequency component
%   OTFo: system OTF
%   co: Wiener filter constant [=1, for minimum RMS estimate]
%   OBJsideP: avg. object spectrum
%   SFo: modulation factor
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
OBJpower = OBJsideP.^2;

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