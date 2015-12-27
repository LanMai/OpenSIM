function OBJparaA = OBJpowerPara(fDIoTnoisy,OTFo)
% AIM: determination of object power parameters Aobj and Bobj
% INPUT VARIABLES
%   fDIoTnoisy: FT of central frequency component
%   OTFo: system OTFo
% OUTPUT VARIABLES
%   OBJparaA: [Aobj Bobj]

w = size(OTFo,1);
wo = w/2;
x = linspace(0,w-1,w);
y = linspace(0,w-1,w);
[X,Y] = meshgrid(x,y);
Cv = (X-wo) + 1i*(Y-wo);
Ro = abs(Cv);

% OTF cut-off frequency
Kotf = OTFedgeF(OTFo);

%% object power parameters through optimization
OBJparaOpt0 = @(OBJpara0)OBJparaOpt(OBJpara0,fDIoTnoisy,OTFo);
options = optimset('LargeScale','off','Algorithm',...
	'active-set','MaxFunEvals',500,'MaxIter',500,'Display','notify');

% obtaining crude initial guesses for Aobj and Bobj 
Zm = (Ro>0.3*Kotf).*(Ro<0.4*Kotf);
Aobj = sum(sum(abs(fDIoTnoisy.*Zm)))./sum(sum(Zm));
Bobj = -0.5;
OBJpara0 = [Aobj Bobj];

% optimization step
[OBJparaA,fval] = fminsearch(OBJparaOpt0,OBJpara0,options);

%% plot for cross-checking the result
%{
Aobj = OBJparaA(1);
Bobj = OBJparaA(2);
OBJo = Aobj*(Ro.^Bobj);
SIGp = OBJo.*OTFo;
pp = 3;
figure;
hold on
plot([0:w-1]-wo,log(abs(fDIoTnoisy(wo+pp,:))),'k--','LineWidth',3,'MarkerSize',6)
plot([0:w-1]-wo,log(abs(SIGp(wo+pp,:))),'mo-','LineWidth',2,'MarkerSize',6)
legend('acutal signal power','avg. signal power')
grid on
box on
%}