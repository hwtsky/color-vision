function [W, paramGA, WI, SX, WP] = Wpig(FC, F, paramGA, d)
% 06/30/2018
% Govardovskii, V. I., et al. In search of the visual pigment template. Visual neuroscience 17, 509-528 (2000).
% for A1 and A2 pigments

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

fmid = 480;
[W1, paramGA, WI1, SX, WP1] = Wpig1(FC, F, paramGA, d);
[W2, paramGA, WI2, SX, WP2] = Wpig2(FC, F, paramGA, d);

indf1 = find(FC>fmid);
indf2 = find(FC<=fmid);

W = [W2(:,indf2), W1(:,indf1)];
WI = [WI2(:,indf2), WI1(:,indf1)];
WP = [WP2(:,indf2), WP1(:,indf1)];

end
