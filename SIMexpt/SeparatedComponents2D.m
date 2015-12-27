function [FiSMao,FiSMap,FiSMam] = SeparatedComponents2D(...
    phaseShift,phaseShift0,FcS1aT,FcS2aT,FcS3aT)
% Aim: Unmixing the frequency components of raw SIM images
%   phaseShift,phaseShift0: illumination phase shifts
%   FcS1aT,FcS2aT,FcS3aT: FT of raw SIM images
%   FiSMao,FiSMap,FiSMam: unmixed frequency components of raw SIM images

phaseShift1 = phaseShift(1,1);
phaseShift2 = phaseShift(2,1);
MF = 1.0;
%% Transformation Matrix
M = 0.5*[1 0.5*MF*exp(-1i*phaseShift0) 0.5*MF*exp(+1i*phaseShift0);
    1 0.5*MF*exp(-1i*phaseShift1) 0.5*MF*exp(+1i*phaseShift1);
    1 0.5*MF*exp(-1i*phaseShift2) 0.5*MF*exp(+1i*phaseShift2)];

%% Separting the components
%===========================================================
Minv = inv(M);

FiSMao = Minv(1,1)*FcS1aT + Minv(1,2)*FcS2aT + Minv(1,3)*FcS3aT;
FiSMap = Minv(2,1)*FcS1aT + Minv(2,2)*FcS2aT + Minv(2,3)*FcS3aT;
FiSMam = Minv(3,1)*FcS1aT + Minv(3,2)*FcS2aT + Minv(3,3)*FcS3aT;
