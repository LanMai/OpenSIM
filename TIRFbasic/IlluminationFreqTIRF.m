function k2fa1 = IlluminationFreqTIRF(fDo,fDp,OTFo,k2o,thetaA)

% AIM: To estimate illumination frequency vector
% INPUT VARIABLES:
%   fDo: noisy central frequency component
%   fDp: noisy off-center frequency component
%   OTFo: system OTF
%   k2o: tentative magnitude of illumination frequency vector
%   thetaA: angular orientation of structured illumination

%% for noise suppression
fAp0 = fDp.*OTFo;
fAo0 = fDo.*OTFo;

% OTF cut-off freq
Kotf = OTFedgeF(OTFo);

w = size(OTFo,1);
wo = w/2;
x = linspace(0,w-1,w);
y = linspace(0,w-1,w);
[X,Y] = meshgrid(x,y);
Cv = (X-wo) + 1i*(Y-wo);
Ro = abs(Cv);
Zo = double( Ro<Kotf );

Ka = k2o.*[cos(thetaA) sin(thetaA)];
Ck = Ka(1) + 1i*Ka(2);
rp = abs(Cv - Ck);
rm = abs(Cv + Ck);
r1 = 20; % local search radius
Zmask = (rp>r1).*(rm>r1);
Zmask = 1 - Zmask;

% computing cross-correlation between central and off-central frequency
% components to find nearest pixel approximation to illumination frequency
% vector
Kmax = k2o + r1;
RHO = zeros(w,w);
for p = wo-Kmax:wo+Kmax
    for q = wo-Kmax:wo+Kmax
        if ( Zmask(p,q) > 0 )
            Cp = (p-wo) + 1i*(q-wo);
            Rp = abs(Cv-Cp);
            Zp = double( Rp<Kotf );
            Z1 = Zo.*Zp;
            fAp1 = circshift(fAp0,[(p-wo) (q-wo)]);
            CC = sum(sum( fAo0.*conj(fAp1) ));
            RHO(p,q) = CC./sum(sum(Z1));
        end
    end
end
% figure, mesh(abs(RHO))

[Mx Idx] = max(abs(RHO),[],1);
[Mx0 Idx0] = max(Mx);
[My Idy] = max(abs(RHO),[],2);
[My0 Idy0] = max(My);
px0 = Idx0 - wo;
py0 = Idy0 - wo;
k2fa = [py0 px0]; % nearest pixel approximation 

%% subpixel approximation
Ifreq2opt0 = @(k2fa0)Ifreq2opt(k2fa0,fDo,fDp,OTFo);
% options = optimset('LargeScale','off','Algorithm',...
% 	'active-set','MaxFunEvals',500,'MaxIter',500,'Display','notify');
options = optimset('LargeScale','off','Algorithm',...
            'active-set','MaxIter',200,'Display','iter');
k2fa0 = k2fa;
[k2fa1,fval] = fminsearch(Ifreq2opt0,k2fa0,options);
k2a = sqrt(k2fa1*k2fa1')
