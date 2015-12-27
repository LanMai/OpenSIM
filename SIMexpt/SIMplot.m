function [] = SIMplot(Fsum,Fperi,Fcent,OTFo,kA,kB,kC,S1aTnoisy)

%% recontructed SIM images
Dcent = real( ifft2(fftshift(Fcent)) );
Dsum = real( ifft2(fftshift(Fsum)) );
Dperi = real( ifft2(fftshift(Fperi)) );

t = size(OTFo,1);
h = 1*30;
figure;
imshow(Dcent(h+1:t-h,h+1:t-h),[])
title('Wiener-Filtered wide-field')
figure;
imshow(Dsum(h+1:t-h,h+1:t-h),[])
title('SIM image')
figure;
imshow(Dperi(h+1:t-h,h+1:t-h),[])
title('SIM image (using only off-center frequency components)')

% appodizing the merged frequency components
Index = 0.4;
Kotf = OTFedgeF(OTFo);
[FsumA] = ApodizationFunction(Fsum,kA,kB,kC,Kotf,Index);
[FperiA] = ApodizationFunction(Fperi,kA,kB,kC,Kotf,Index);
DsumA = real( ifft2(fftshift(FsumA)) );
DperiA = real( ifft2(fftshift(FperiA)) );

figure;
imshow(DsumA(h+1:t-h,h+1:t-h),[])
title('appodized SIM image')
figure;
imshow(DperiA(h+1:t-h,h+1:t-h),[])
title('appodized SIM image (using only off-center frequency components)')


%% Frequency Plots
fS1aTnoisy = fftshift(fft2(S1aTnoisy));
w = size(OTFo,1);
wo = w/2;
w1 = size(fS1aTnoisy,1);
if ( w > w1 )
    fS1aTnoisy0 = zeros(w,w);
    fS1aTnoisy0(wo-w1/2+1:wo+w1/2,wo-w1/2+1:wo+w1/2) = fS1aTnoisy;
    clear fS1aTnoisy
    fS1aTnoisy = fS1aTnoisy0;
    clear fS1aTnoisy0;
end


p = 10;
minL1 = min(min( abs(fS1aTnoisy).^(1/p) ));
minL2 = min(min( abs(Fcent).^(1/p) ));
minL3 = min(min( abs(Fsum).^(1/p) ));
minL4 = min(min( abs(Fperi).^(1/p) ));
minL5 = min(min( abs(FsumA).^(1/p) ));
maxL1 = max(max( abs(fS1aTnoisy).^(1/p) ));
maxL2 = max(max( abs(Fcent).^(1/p) ));
maxL3 = max(max( abs(Fsum).^(1/p) ));
maxL4 = max(max( abs(Fperi).^(1/p) ));
maxL5 = max(max( abs(FsumA).^(1/p) ));
minL = min([minL1,minL2,minL3,minL4,minL5]);
maxL = max([maxL1,maxL2,maxL3,maxL4,maxL5]);

figure;
imshow(abs(fS1aTnoisy).^(1/p),[minL maxL])
title('fS1aTnoisy')

figure;
imshow(abs(Fcent).^(1/p),[minL maxL])
title('Weiner Filtered frequency')
figure;
imshow(abs(Fsum).^(1/p),[minL maxL])
title('SIM frequency')
figure;
imshow(abs(Fperi).^(1/p),[minL maxL])
title('SIM (off-center frequency components)')
figure;
imshow(abs(FsumA).^(1/p),[minL maxL])
title('appodized SIM frequency')
figure;
imshow(abs(FperiA).^(1/p),[minL maxL])
title('appodized SIM (off-center frequency components)')


