function [Fsum,Fperi,Fcent] = MergingHeptaletsF(fDIo,fDIp,fDIm,fDBo,fDBp,fDBm,...
    fDCo,fDCp,fDCm,Ma,Mb,Mc,npDIo,npDIp,npDIm,...
    npDBo,npDBp,npDBm,npDCo,npDCp,npDCm,k2fa,k2fb,k2fc,OBJpara,OTFo)

% AIM: To merge all 9 frequency components into one using generalized
%       Wiener Filter
% INPUT VARIABLES
%   [fDIo fDIp fDIm
%    fDBo fDBp fDBm
%    fDCo fDCp fDCm]: nine frequency components
%   Ma,Mb,Mc: modulation factors for the three illumination orientations
%   [npDIo npDIp npDIm 
%    npDBo npDBp npDBm 
%    npDCo,npDCp,npDCm]: noise powers corresponding to nine frequency
%                       components
%   k2fa,k2fb,k2fc: illumination frequency vectors for the three
%                   illumination orientations
%   OBJpara: Object spectrum parameters
%   OTFo: system OTF
% OUTPUT VARIABLES
%   Fsum: all nine frequency components merged into one using 
%           generalised Wiener Filter
%   Fperi: six off-center frequency components merged into one using 
%           generalised Wiener Filter
%   Fcent: averaged of the three central frequency components


% obtain signal spectrums corresponding to central and off-center
% frequency components
[SIGao,SIGam2,SIGap2] = TripletSNR0(OBJpara,k2fa,OTFo,fDIm);
[SIGbo,SIGbm2,SIGbp2] = TripletSNR0(OBJpara,k2fb,OTFo,fDBm);
[SIGco,SIGcm2,SIGcp2] = TripletSNR0(OBJpara,k2fc,OTFo,fDCm);

SIGap2 = Ma*SIGap2;
SIGam2 = Ma*SIGam2;
SIGbp2 = Mb*SIGbp2;
SIGbm2 = Mb*SIGbm2;
SIGcp2 = Mc*SIGcp2;
SIGcm2 = Mc*SIGcm2;

%% Generalized Wiener-Filter computation
SNRao = SIGao.*conj(SIGao)./npDIo;
SNRap = SIGap2.*conj(SIGap2)./npDIp;
SNRam = SIGam2.*conj(SIGam2)./npDIm;

SNRbo = SIGbo.*conj(SIGbo)./npDBo;
SNRbp = SIGbp2.*conj(SIGbp2)./npDBp;
SNRbm = SIGbm2.*conj(SIGbm2)./npDBm;

SNRco = SIGco.*conj(SIGco)./npDCo;
SNRcp = SIGcp2.*conj(SIGcp2)./npDCp;
SNRcm = SIGcm2.*conj(SIGcm2)./npDCm;

ComDeno = 0.01 + ( SNRao + SNRap + SNRam + SNRbo + SNRbp + SNRbm + SNRco + SNRcp + SNRcm );
Fsum = fDIo.*SNRao + fDIp.*SNRap + fDIm.*SNRam...
    + fDBo.*SNRbo + fDBp.*SNRbp + fDBm.*SNRbm...
    + fDCo.*SNRco + fDCp.*SNRcp + fDCm.*SNRcm;
Fsum = Fsum./ComDeno;

ComPeri = 0.01 + ( SNRap + SNRam + SNRbp + SNRbm + SNRcp + SNRcm );
Fperi = fDIp.*SNRap + fDIm.*SNRam + fDBp.*SNRbp + fDBm.*SNRbm...
    + fDCp.*SNRcp + fDCm.*SNRcm;
Fperi = Fperi./ComPeri;

% averaged central frequency component
Fcent = (fDIo+fDBo+fDCo)/3;


