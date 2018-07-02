function [W, paramGA, WI, SX, WP] = Wpig1(FC, F, paramGA, d)
% 06/30/2018
% Govardovskii, V. I., et al. In search of the visual pigment template. Visual neuroscience 17, 509-528 (2000).
% for A1 pigments

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

b = 0.922;
c = 1.104;
A = 69.7;
B = 28;
C = -14.9;
D = 0.674;

Ea = 0.8795 + 0.0459*exp(-((FC-300).^2)/11940);
FD = repmat(FC', KF, 1)./repmat(F, 1, NA);

W = 1./(exp(A*(repmat(Ea',KF,1)-FD))+ exp(B*(b-FD)) + exp(C*(c-FD))+D);

Lb = 189 + 0.315*FC;
Eb = -40.5 + 0.195*FC;

Sb = 0.26*exp(-((repmat(F, 1, NA)-repmat(Lb',KF,1))./repmat(Eb',KF,1)).^2);

W = W + Sb;

indx = 1:KF1;%

WI = W(len0+1:len0+KF1,:)./repmat(max(W(len0+1:len0+KF1,:)), KF, 1);
%-----------------------------------
WP = 0;
if paramGA.Dpeak
    W = 1 - 10.^(-paramGA.Dpeak*W);
    WP = W(len0+1:len0+KF1,:);
    if ~paramGA.Trans
        if paramGA.FgN
            W = W./repmat(max(W), KF, 1);
        end
    end
end

SX = FitLenMacu(F, Lc, Mc, paramGA.isPoly);
if paramGA.Trans
    W = W.*repmat(SX, 1,NA);
    if paramGA.FgN
        W = W./repmat(max(W), KF, 1);
    end
end
SX = SX(len0+1:len0+KF1,:);
SX = SX(indx,:);

end





