function [W, paramGA, WI, SX, WP] = Wpig2(FC, F, paramGA, d)
% 06/30/2018
% Govardovskii, V. I., et al. In search of the visual pigment template. Visual neuroscience 17, 509-528 (2000).
% for A2 pigments

if ~exist('FC','var')
    FC = [350:1:565]';
end

if ~exist('F','var')
    F = [380:1:700]';
end

if ~exist('d','var')
    d = round(F(2) - F(1));
end

if ~exist('paramGA','var')
    paramGA.r = 0;
    paramGA.Dpeak = 0.3;
    paramGA.Trans = 1;
    paramGA.FgN = 0;
    paramGA.isPoly = 0;
    paramGA.Lc = 1.0;
    paramGA.Mc = 1.0;
end

Lc = paramGA.Lc;
Mc = paramGA.Mc;

NA = length(FC);

Fmin = min(F);
Fmax = max(F);
df0 = F(2)-F(1);
df = df0/d;
Fd = F;%[Fmin:df:Fmax]';
KF1 = length(Fd);

r = paramGA.r;
F0 = [Fmin-r*50:df:Fmin-df]';
F1 = [Fmax+df:df:Fmax+r*50]';

F = [F0;Fd;F1];

len0 = length(F0);
KF = length(F);

b = 0.9101;
c = 1.1123;
B = 20.85;
C = -10.37;
D = 0.5343;

Ea = 0.875 + 0.0268*exp(((FC-665))/40.7);%a;%
EA = 62.7 + 1.834*exp((FC-625)/54.2);%A;%
FD = repmat(FC', KF, 1)./repmat(F, 1, NA);%

W = 1./(exp(repmat(EA',KF,1).*(repmat(Ea',KF,1)-FD))+ exp(B*(b-FD)) + exp(C*(c-FD))+D) ;%

Lb = 216.7 + 0.287*FC;
Eb = 317 - 1.149*FC + 0.00124*FC.*FC;

Sb = 0.37*exp(-((repmat(F, 1, NA) - repmat(Lb',KF,1))./repmat(Eb',KF,1)).^2);

W = W + Sb;

indx = 1:KF1;

WI = W(len0+1:len0+KF1,:)./repmat(max(W(len0+1:len0+KF1,:)), KF, 1);%;%
%-----------------------------------
WP = 0;
if paramGA.Dpeak
    W = 1 - 10.^(-paramGA.Dpeak*W);
    WP = W(len0+1:len0+KF1,:);
    if ~paramGA.Trans
        if paramGA.FgN
            W = W./repmat(max(W), KF, 1);%
        end
    end
end

SX = FitLenMacu(F, Lc, Mc, paramGA.isPoly);
if paramGA.Trans
    W = W.*repmat(SX, 1,NA);
    if paramGA.FgN
        W = W./repmat(max(W), KF, 1);%
    end
end
SX = SX(len0+1:len0+KF1,:);
SX = SX(indx,:);

end

