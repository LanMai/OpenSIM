function CCo = Kai2Opt(phaseShift,FcS1aT,FcS2aT,OTFo)
% Aim: Computes normalized cross-power spectrum between separated frequency
% components
% INPUT VARIABLES: 
%   phaseShift: relative phases [phaseShift2,phaseShift3]
%   FcS1aT: FT of difference between raw images 1 and 2
%   FcS2aT: FT of difference between raw images 1 and 3
%   OTFo: system OTF

phaseShift2 = phaseShift(1);
phaseShift3 = phaseShift(2);

phaseShift1 = 0;
phase2 = exp(-1i*phaseShift2) - exp(-1i*phaseShift1);
phase3 = exp(-1i*phaseShift3) - exp(-1i*phaseShift1);
%% Transformation Matrix
M = [phase2 conj(phase2);
    phase3 conj(phase3)];

%% Separting the components
%===========================================================
Minv = inv(M);

FiSMap = Minv(1,1)*FcS1aT + Minv(1,2)*FcS2aT;
FiSMam = Minv(2,1)*FcS1aT + Minv(2,2)*FcS2aT;

%{
figure;
mesh(log(abs(FiSMap)))
figure;
mesh(log(abs(FiSMam)))
kk
%}

FiSMap0 = FiSMap.*conj(OTFo);
FiSMam0 = FiSMam.*conj(OTFo);

CCo = abs( sum(sum( FiSMap0.*conj(FiSMam0) )) );
