function [S1aTnoisyA, S2aTnoisyA, S3aTnoisyA, ...
 S1bTnoisyA, S2bTnoisyA, S3bTnoisyA, ...
 S1cTnoisyA, S2cTnoisyA, S3cTnoisyA] = BgNormalization0(...
    S1aTnoisy, S2aTnoisy, S3aTnoisy, ...
    S1bTnoisy, S2bTnoisy, S3bTnoisy, ...
    S1cTnoisy, S2cTnoisy, S3cTnoisy) 

m1a = mean2(S1aTnoisy);
m2a = mean2(S2aTnoisy);
m3a = mean2(S3aTnoisy);
m1b = mean2(S1bTnoisy);
m2b = mean2(S2bTnoisy);
m3b = mean2(S3bTnoisy);
m1c = mean2(S1cTnoisy);
m2c = mean2(S2cTnoisy);
m3c = mean2(S3cTnoisy);

s1a = std(S1aTnoisy(:));
s2a = std(S2aTnoisy(:));
s3a = std(S3aTnoisy(:));
s1b = std(S1bTnoisy(:));
s2b = std(S2bTnoisy(:));
s3b = std(S3bTnoisy(:));
s1c = std(S1cTnoisy(:));
s2c = std(S2cTnoisy(:));
s3c = std(S3cTnoisy(:));

mTemp = max([m1a m2a m3a m1b m2b m3b m1c m2c m3c]);
sTemp = max([s1a s2a s3a s1b s2b s3b s1c s2c s3c]);

S1aTnoisyA = (S1aTnoisy - m1a).*(sTemp./s1a) + mTemp;
S2aTnoisyA = (S2aTnoisy - m2a).*(sTemp./s2a) + mTemp;
S3aTnoisyA = (S3aTnoisy - m3a).*(sTemp./s3a) + mTemp;
S1bTnoisyA = (S1bTnoisy - m1b).*(sTemp./s1b) + mTemp;
S2bTnoisyA = (S2bTnoisy - m2b).*(sTemp./s2b) + mTemp;
S3bTnoisyA = (S3bTnoisy - m3b).*(sTemp./s3b) + mTemp;
S1cTnoisyA = (S1cTnoisy - m1c).*(sTemp./s1c) + mTemp;
S2cTnoisyA = (S2cTnoisy - m2c).*(sTemp./s2c) + mTemp;
S3cTnoisyA = (S3cTnoisy - m3c).*(sTemp./s3c) + mTemp;