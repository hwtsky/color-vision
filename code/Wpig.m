function [W, paramGA, WI, SX, WP] = Wpig(FC, F, paramGA, d)
% 06/30/2018
% Govardovskii, V. I., et al. In search of the visual pigment template. Visual neuroscience 17, 509-528 (2000).
% for A1 and A2 pigments

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

fmid = 480;
[W1, paramGA, WI1, SX, WP1] = Wpig1(FC, F, paramGA, d);
[W2, paramGA, WI2, SX, WP2] = Wpig2(FC, F, paramGA, d);

indf1 = find(FC>fmid);
indf2 = find(FC<=fmid);

W = [W2(:,indf2), W1(:,indf1)];
WI = [WI2(:,indf2), WI1(:,indf1)];
WP = [WP2(:,indf2), WP1(:,indf1)];

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

