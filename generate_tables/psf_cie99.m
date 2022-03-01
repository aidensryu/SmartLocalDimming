function psf=psf_cie99(angledeg, age, p)
% psf=psf_cie99(angledeg, age, p)
%
% angledeg - angle in visual degrees
% age - age in years
% p - iris pigmentation
%  p=0   - Non caucasian
%  p=0.5 - Brown eye
%  p=1   - Green eye
%  p=1.2 - Light blue
%
% This model does not account for a pupil size (assumed that it is readly
% known).
%
% This is a glare spread function from:
% Vos, J.J. and van den Berg, T.J.T.P [CIE Research note 135/1, “Disability
% Glare”] ISBN 3 900 734 97 6 (1999).

if nargin<2
    age=0;
end

if nargin<3
    %p=0     % Non caucasian
    %p=0.5   % Brown eye
    p=1;     % Green eye
    %p=1.2   % Light blue
end

R=angledeg;
agefac=(age/70)^4;
% core
psf=(1-0.08*agefac)*(9.2e6./(1+(R/0.0046).^2).^1.5 + 1.5e5./(1+(R/0.045).^2).^1.5);
% skirt
psf=psf+(1+1.6*agefac)*(400./(1+(R/0.1).^2) + 3e-8*R.^2);
% PE
psf=psf+p*(1+1.6*agefac)*(1300./(1+(R/0.1).^2).^1.5 + 0.8./(1+(R/0.1).^2).^0.5);
% wall
psf=psf+p*2.5e-3;

normcorr=p*(0.0417+agefac*0.055);

psf=psf/(1+normcorr);

