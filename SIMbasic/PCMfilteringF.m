function [fDof,fDp2,fDm2,npDo,npDp,npDm,Mm,DoubleMatSize]...
    = PCMfilteringF(fDo,fDp,fDm,OTFo,OBJparaA,kA)

% AIM: obtaining Wiener Filtered estimates of noisy frequency components
% INPUT VARIABLES
%   fDo,fDp,fDm: noisy estimates of separated frequency components
%   OTFo: system OTF
%   OBJparaA: object power parameters
%   kA: illumination vector
% OUTPUT VARIABLES
%   fDof,fDp2,fDm2: Wiener Filtered estimates of fDo,fDp,fDm
%   npDo,npDp,npDm: avg. noise power in fDo,fDp,fDm
%   Mm: illumination modulation factor
%   DoubleMatSize: parameter for doubling FT size if necessary

w = size(OTFo,1);
wo = w/2;
x = linspace(0,w-1,w);
y = linspace(0,w-1,w);
[X,Y] = meshgrid(x,y);
Cv = (X-wo) + 1i*(Y-wo);
Ro = abs(Cv);

% OTF cut-off frequency
Kotf = OTFedgeF(OTFo);

% object power parameters
Aobj = OBJparaA(1);
Bobj = OBJparaA(2);

% Wiener Filtering central frequency component
SFo = 1;
co = 1.0; 
[fDof,npDo] = WoFilterCenter(fDo,OTFo,co,OBJparaA,SFo);

% modulation factor determination
Mm = ModulationFactor(fDp,kA,OBJparaA,OTFo);

%% Duplex power (default)
kv = kA(2) + 1i*kA(1); % vector along illumination direction
Rp = abs(Cv-kv);
Rm = abs(Cv+kv);
OBJp = Aobj*(Rp.^Bobj);
OBJm = Aobj*(Rm.^Bobj);
k3 = round(kA);
OBJp(wo+1+k3(1),wo+1+k3(2)) = 0.25*OBJp(wo+2+k3(1),wo+1+k3(2))...
	+ 0.25*OBJp(wo+1+k3(1),wo+2+k3(2))...
	+ 0.25*OBJp(wo+0+k3(1),wo+1+k3(2))...
	+ 0.25*OBJp(wo+1+k3(1),wo+0+k3(2));
OBJm(wo+1-k3(1),wo+1-k3(2)) = 0.25*OBJm(wo+2-k3(1),wo+1-k3(2))...
	+ 0.25*OBJm(wo+1-k3(1),wo+2-k3(2))...
	+ 0.25*OBJm(wo+0-k3(1),wo+1-k3(2))...
	+ 0.25*OBJm(wo+1-k3(1),wo+0-k3(2));

% Filtering side lobes (off-center frequency components)
SFo = Mm;
[fDpf,npDp] = WoFilterSideLobe(fDp,OTFo,co,OBJm,SFo);
[fDmf,npDm] = WoFilterSideLobe(fDm,OTFo,co,OBJp,SFo);

%% doubling Fourier domain size if necessary
DoubleMatSize = 0;
if ( 2*Kotf > wo )
	DoubleMatSize = 1; % 1 for doubling fourier domain size, 0 for keeping it unchanged
end
if ( DoubleMatSize>0 )
    t = 2*w;
    to = t/2;
    u = linspace(0,t-1,t);
    v = linspace(0,t-1,t);
    [U,V] = meshgrid(u,v);
    fDoTemp = zeros(2*w,2*w);
    fDpTemp = zeros(2*w,2*w);
    fDmTemp = zeros(2*w,2*w);
    OTFtemp = zeros(2*w,2*w);
    fDoTemp(wo+1:w+wo,wo+1:w+wo) = fDof;
    fDpTemp(wo+1:w+wo,wo+1:w+wo) = fDpf;
    fDmTemp(wo+1:w+wo,wo+1:w+wo) = fDmf;
    OTFtemp(wo+1:w+wo,wo+1:w+wo) = OTFo;
    clear fDof fDpf fDmf OTFo w wo x y X Y
    fDof = fDoTemp;
    fDpf = fDpTemp;
    fDmf = fDmTemp;
    OTFo = OTFtemp;
    clear fDoTemp fDpTemp fDmTemp OTFtemp
else
    t = w;
    to = t/2;
    u = linspace(0,t-1,t);
    v = linspace(0,t-1,t);
    [U,V] = meshgrid(u,v);
end

% Shifting the off-center frequency components to their correct location
fDp1 = fft2(ifft2(fDpf).*exp( +1i.*2*pi*(kA(2)/t.*(U-to) + kA(1)/t.*(V-to)) ));
fDm1 = fft2(ifft2(fDmf).*exp( -1i.*2*pi*(kA(2)/t.*(U-to) + kA(1)/t.*(V-to)) ));

%% Shift induced phase error correction
Cv = (U-to) + 1i*(V-to);
Ro = abs(Cv);
Rp = abs(Cv-kv);
k2 = sqrt(kA*kA');

% frequency range over which corrective phase is determined
Zmask = (Ro < 0.8*k2).*(Rp < 0.8*k2);

% corrective phase
Angle0 = angle( sum(sum( fDof.*conj(fDp1).*Zmask )) );

% phase correction
fDp2 = exp(+1i*Angle0).*fDp1;
fDm2 = exp(-1i*Angle0).*fDm1;

%% for visual verification
%{
pp = 3;
figure;
hold on
plot(u-to,angle(fDof(to+pp,:)).*180/pi,'o-','LineWidth',2,'MarkerSize',8)
plot(u-to,angle(fDp1(to+pp,:)).*180/pi,'ro-','LineWidth',2,'MarkerSize',8)
plot(u-to,angle(fDp2(to+pp,:)).*180/pi,'go--','LineWidth',2,'MarkerSize',6)
grid on
box on
figure;
hold on
plot3(u-to,real(fDof(to+1,:)),imag(fDof(to+1,:)),'o-','LineWidth',2,'MarkerSize',8)
plot3(u-to,real(fDp1(to+1,:)),imag(fDp1(to+1,:)),'ro-','LineWidth',2,'MarkerSize',8)
plot3(u-to,real(fDp2(to+1,:)),imag(fDp2(to+1,:)),'go--','LineWidth',2,'MarkerSize',6)
grid on
box on

ff1 = fDof.*conj(fDp1);
ff2 = fDof.*conj(fDp2);
figure;
hold on
plot(u-to,angle(ff1(to+pp,:)).*180/pi,'ro-','LineWidth',2,'MarkerSize',8)
plot(u-to,angle(ff2(to+pp,:)).*180/pi,'go--','LineWidth',2,'MarkerSize',6)
grid on
box on
%}




