function Kotf = OTFedgeF(OTFo)

w = size(OTFo,1);
wo = w/2;

OTF1 = OTFo(wo+1,:);
OTFmax = max(max(abs(OTFo)));
OTFtruncate = 0.01;
i = 1;
while ( abs(OTF1(1,i))<OTFtruncate*OTFmax )
	Kotf = wo+1-i;
	i = i + 1;
end 