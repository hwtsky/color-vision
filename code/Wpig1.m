function [W, paramGA, WI, SX, WP] = Wpig1(FC, F, paramGA, d)
% Govardovskii, V. I., et al. In search of the visual pigment template. Visual neuroscience 17, 509-528 (2000).
% for A1 pigments

if ~exist('FC','var')    
    FC = [380:1:700]';%[410:700]';
end

if ~exist('F','var')
    F = [380:1:700]';%[410:10:700]';
end

if ~exist('d','var')
    d = round(F(2) - F(1));
end

% LamF = sqrt(F'*F);
% FL = F/LamF;

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

b = 0.922;
c = 1.104;
A = 69.7;
B = 28;
C = -14.9;
D = 0.674;

Ea = 0.8795 + 0.0459*exp(-((FC-300).^2)/11940);%a;%
FD = repmat(FC', KF, 1)./repmat(F, 1, NA);%

W = 1./(exp(A*(repmat(Ea',KF,1)-FD))+ exp(B*(b-FD)) + exp(C*(c-FD))+D) ;%

Lb = 189 + 0.315*FC;
Eb = -40.5 + 0.195*FC;

Sb = 0.26*exp(-((repmat(F, 1, NA)-repmat(Lb',KF,1))./repmat(Eb',KF,1)).^2);

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
%     filename = ['../data/TransFun.mat'];
%     load(filename, 'Lambda', 'XL', 'XM', 'PLens', 'PMacular', 'PLenMac', 'LenCut');
%     SLM = polyval(PLenMac, F);%PLenMac 

    W = W.*repmat(SX, 1,NA);
    
%     b = 0.01;
%     a = 400;
%     T = 1-(exp(-b*(FC-a)));%exp-
%     W = W.*repmat(T', KF,1);  
    
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



