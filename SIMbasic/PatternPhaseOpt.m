function CCop = PatternPhaseOpt(phaseA,S1aTnoisy,k2fa)

w = size(S1aTnoisy,1);
wo = w/2;

x = linspace(0,w-1,w);
y = linspace(0,w-1,w);
[X,Y] = meshgrid(x,y);

sAo = cos( 2*pi*(k2fa(2).*(X-wo)+k2fa(1).*(Y-wo))./w + phaseA );
S1aTnoisy = S1aTnoisy - mean2(S1aTnoisy);
CCop = -sum(sum(S1aTnoisy.*sAo));