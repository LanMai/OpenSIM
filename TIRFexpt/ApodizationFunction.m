function [FsumA] = ApodizationFunction(Fsum,k2fa,k2fb,k2fc,Kotf,Index)

t = size(Fsum,1);

OTFo_mask = OTFmaskShifted(0.*k2fa,Kotf,t);
OTFap_mask = OTFmaskShifted(k2fa,Kotf,t);
OTFam_mask = OTFmaskShifted(-k2fa,Kotf,t);
OTFbp_mask = OTFmaskShifted(k2fb,Kotf,t);
OTFbm_mask = OTFmaskShifted(-k2fb,Kotf,t);
OTFcp_mask = OTFmaskShifted(k2fc,Kotf,t);
OTFcm_mask = OTFmaskShifted(-k2fc,Kotf,t);
ApoMask = OTFo_mask.*OTFap_mask.*OTFam_mask.*OTFbp_mask.*OTFbm_mask.*OTFcp_mask.*OTFcm_mask;
clear OTFoz OTFo_mask OTFap_mask OTFam_mask OTFbp_mask OTFbm_mask OTFcp_mask OTFcm_mask
%{
figure;
imshow(ApoMask,[ ])
%}

DistApoMask = bwdist(ApoMask);
maxApoMask = max(max(DistApoMask));
ApoFunc = double(DistApoMask./maxApoMask).^Index;
%{
figure;
imshow(DistApoMask,[ ])
figure;
mesh(double(ApoFunc))
%}

FsumA = Fsum.*ApoFunc;
clear ApoMask DistApoMask ApoFunc