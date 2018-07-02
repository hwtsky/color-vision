function [W, paramGA, WI, SX, WP] = Wpig2(FC, F, paramGA, d)
% 06/30/2018
% Govardovskii, V. I., et al. In search of the visual pigment template. Visual neuroscience 17, 509-528 (2000).
% for A2 pigments

if ~exist('FC','var')
    FC = [380:1:700]';%[410:700]';
end

if ~exist('F','var')
    F = [380:1:700]';%[410:10:700]';
end

if ~exist('d','var')
    d = round(F(2) - F(1));
end

if ~exist('paramGA','var')
    paramGA.r = 0;
    paramGA.Dpeak = 0.3;%0.3
    paramGA.Trans = 1;
    paramGA.FgN = 0;
    paramGA.isPoly = 0;
    paramGA.Lc = 1.0;
    paramGA.Mc = 1.0;
end

Lc = paramGA.Lc;
Mc = paramGA.Mc;

NA = length(FC);
KF0 = length(F);

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

W = W + Sb;%

indx = 1:KF1;%d:

W0 = W;

%-----------------------------------

WI = W(len0+1:len0+KF1,:)./repmat(max(W(len0+1:len0+KF1,:)), KF, 1);%;%
%-----------------------------------
WP = 0;
if paramGA.Dpeak
    W = 1 - 10.^(-paramGA.Dpeak*W);
    WP = W(len0+1:len0+KF1,:);
    if ~paramGA.Trans
        if paramGA.FgN
            W = W./repmat(max(W), KF, 1);%sum(W) %max(W)
        end
    end
end

SX = FitLenMacu(F, Lc, Mc, paramGA.isPoly);%_2
if paramGA.Trans
    W = W.*repmat(SX, 1,NA);
    if paramGA.FgN
        W = W./repmat(max(W), KF, 1);%sum(W) %max(W) 
    end
end
SX = SX(len0+1:len0+KF1,:);
SX = SX(indx,:);
%-----------------------------------------------
Wd = (W(2:end,:)-W(1:end-1,:))/df;
Wd(end+1,:) = (W(end,:)-W(end-1,:))/df;
W = W(len0+1:len0+KF1,:);
Wd = Wd(len0+1:len0+KF1,:);

W = W(indx,:);
Wd = Wd(indx,:);

return;

%%

if ~exist('rw','var')
    rw = 0.9999;
end

WL = W;%log10(W);
[U,D,Vd] = svd([W]);%
d = diag(D);
d2 = d.*d;
ds = sqrt(cumsum(d2));
dc = ds./sqrt(sum(d2));
ind = find(dc<=rw);%
len = length(ind);

UT = U(:,1:len);
WT = UT*UT'*W;

%---------------------------------

F = F(len0+1:len0+KF1);
F = F(indx);
ID = zeros(NA,1);
ID(1) = 1;
ID(end) = 1;
ID((floor(FC-430) == 0) | (floor(FC-530) == 0)| (floor(FC-560) == 0)) = 1;

% figure, 
% plot(FC,T);

figure,
plot(F,SX)

d = 50;

% figure('Name',['PigTempFB4, W']);
% hold on
% for j = 1:d:NA
%     if 1%ID(j)%
%         w = W0(:,j);
%         plot(F,(w),'b')
%     end
% end


figure('Name',['PigTempFB4, WT: ', num2str(len), '_', num2str(rw)]);
hold on
for j = 1:d:NA
    if 1%ID(j)%
        w = W(:,j);
        plot(F,w,'r')
    end
end

% figure('Name',['PigTempFB4, Wd']);
% hold on
% for j = 1:d:NA
%     if ID(j)%1%
%         w = Wd(:,j);
%         plot(F,(w),'b')
%     end
% end

% figure('Name',['PigTempFB4, Wd*FL']);
% hold on
% for j = 1:d:NA
%     if ID(j)%1%
%         w = Wd(:,j);
%         plot(F,(w.*FL),'r')
%     end
% end

% figure('Name',['W_Log']);
% hold on
% for j = 1:50:NA
%     w = W(:,j);
%     plot(F,log10(w),'b')
% end

% if paramGA.Sig
%     plot(F,B,'r--')
% end

