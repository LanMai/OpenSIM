function OTF2 = OTFdoubling(OTFo,DoubleMatSize)

%% embeds OTFo in doubled range of frequencies  

w = size(OTFo,1);
wo = w/2;
if ( DoubleMatSize>0 )
    t = 2*w;
    u = linspace(0,t-1,t);
    v = linspace(0,t-1,t);
    OTF2 = zeros(2*w,2*w);
    OTF2(wo+1:w+wo,wo+1:w+wo) = OTFo;
    clear OTFo
else
    OTF2 = OTFo;
    clear OTFo
end