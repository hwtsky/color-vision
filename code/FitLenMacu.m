function [SX, FLenMac] = FitLenMacu(F, Lc, Mc, isPoly)
%% FitLenMacu
if ~exist('F', 'var')
    F =[380:1:700]';
end
if ~exist('Lc', 'var')
    Lc = 1.0;
end
if ~exist('Mc', 'var')
    Mc = 1.0;
end
if ~exist('isPoly', 'var')
    isPoly = 0;
end

Lambda = [380:5:800]';
N = length(Lambda);

Lens = zeros(N,1);
Macular = zeros(N,1);

XL = [2.0, 1.8, 1.6, 1.4,...
    1.2, 1.0, 0.82, 0.67, 0.55, 0.45, 0.37, 0.31, 0.27, 0.24, 0.225, 0.205,...
    .195, .185, .175, .165, .155, .145, .14, .13, .125, .115, .11, .105, .1,...
    .095, .09,  .085, .08, .075, .07, .065, .06, .055, .05, .045, .04, .035,...
    .031, .027, .024, .021, .018, .015, .012, .01, .008, .006, .004, .003, .002, .001]';
XM = [0, 0, 0.015, .05,...
    .085, .12, .16, .225, .3, .345, .365, .38, .4, .425, .46, .49, .495, .47, ...
    .445, .41, .415, .42, .41, .36, .275, .195, .13, .085, .05, .025, .01];


Lens(1:length(XL)) = XL;
Macular(1:length(XM)) = XM;

LenCut = Lambda(length(XL));
MacuCut = Lambda(length(XM));

if isPoly == 0
    FLens = interp1(Lambda,Lens,F);%
    FMacular = interp1(Lambda,Macular,F);
    FLenMac = Lc*FLens + Mc*FMacular;
else
    np = 20;
    PLens = polyfit(Lambda,Lens,np);
    PMacular = polyfit(Lambda,Macular,np);
    FLens = polyval(PLens,F);
    FMacular = polyval(PMacular,F);
    FLenMac = Lc*FLens + Mc*FMacular;
end

FLenMac(FLenMac<0) = 0;

SX = 10.^(-FLenMac);%
SX(F>max([LenCut, MacuCut])) = 1;
SX(F<360) = 0;

return;

