%% Main 06/30/2018
%==========================================================================
clc
close all
clear all

%% Init Parameters
fprintf('\nStart Main ...\n');
%--------------------------------------------------------------------------
pdfname = 'CDNS'; % 'CDNS' 'Forest' 'HINS' 'Munsell'
typeNo = 0; %0;% 1;% 2;% 3;%
%------------------------------------- 
fileNo = '1-1';%'1-2';%'2-1';%'2-2';%

postfix = 'c';

a = 170;%100;%130;%%170;%1000;%5000;%10000;

% a = 200;%2000;% 

%-------------------------------------
MaxIter = 1000;
SampleNum = 10000;
epsOP = 1e-300;
seed = typeNo + 1;
%-------------------------------------
KX = 10;%20;%
DF = 1;
DFC = 1;%0.5;%

Fmin = 380;% >= 380
Fmax = 700;

FCmin = 350;%Fmin;% <= Fmin
FCmax = 565;% 560 565 570 580 600
%-------------------------------------
paramGA.Dpeak = 0.3;
paramGA.Trans = 1;
paramGA.FgN = 0;
paramGA.r = 0;
paramGA.Lc = 1.0;%
paramGA.Mc = 1.0;%
paramGA.isPoly = 0;

%--------------------------------------------------------------------------
%% Load data
F = [Fmin:DF:Fmax]';
FC = [FCmin:DFC:FCmax]';

NA = length(FC);

[W, paramGA, WI] = Wpig(FC, F, paramGA);

%--------------------------------------------------------------------------
fn = [pdfname, '_',num2str(typeNo),'_',num2str(fileNo), ...
      '_', num2str(FCmax),'_', num2str(a)];
figurename = ['Fig', num2str(fileNo), '_', postfix];%fn;%

FigName = ['../results/figures/', figurename];

FontSize = 10;
LineWidth = 0.5;

%% --------------------------------------------------------------------------
if ~ismember({'X'},who) || isempty(X)
    pdfparams.typeNo = typeNo;
    pdfparams.Fmin = Fmin;
    pdfparams.Fmax = Fmax;
    pdfparams.seed = seed;
    pdfparams.SampleNum = SampleNum;
    if strcmp(pdfname, 'CDNS') 
        if typeNo == 0
            pdfparams.reflectfn = {'ugmlrf'};
        elseif typeNo == 1
            pdfparams.reflectfn = {'ugylrf'};
        elseif typeNo >= 2
            pdfparams.r1 = 1;
            pdfparams.r2 = 1;
            pdfparams.reflectfn = {'ugmlrf', 'ugylrf'};
        end
    end
    
    NX = length(F);
	[X, F0] = GetDataset(pdfname, pdfparams);
    [NX0, SampleNum] = size(X);
    if (NX0 ~= NX) || (Fmin ~= F0(1)) || (Fmax ~= F0(end))
        X0 = zeros(NX, SampleNum);
        for n = 1:SampleNum
            f = griddedInterpolant(F0, X(:,n),'cubic');%'linear'
            X0(:,n) = f(F);
        end
        X = X0;
        X0 = [];
    end
    X(X<0) = 0;
    X = X/sum(mean(X,2));
end
%% 
[NX, SampleNum] = size(X);

Param.pdfname = pdfname;
Param.typeNo = typeNo;
Param.fileNo = fileNo;
Param.a = a;
Param.NA = NA;
Param.NX = NX;
Param.KX = KX;
Param.MaxIter = MaxIter;
Param.SampleNum = SampleNum;
Param.Fmin = Fmin;
Param.Fmax = Fmax;
Param.FCmin = FCmin;
Param.FCmax = FCmax;
Param.DF = DF;
Param.DFC = DFC;
Param.seed = seed;
Param.epsOP = epsOP;
Param.paramGA = paramGA;
Param.pdfparams = pdfparams;
%--------------------------------------------------------------------------
%% SVD

xv = mean(X, 2);
X = X - repmat(xv, 1, SampleNum);
[U, D, Vd] = svd(X);
Vd = [];
d = diag(D);

% ---------------------------
Wu = U'*W;
sigx = d/sqrt(SampleNum);
DSigx = diag(sigx);
Ksx = length(sigx);
WN = a*DSigx*Wu(1:Ksx,:);
[Uw,Dw,Vdw] = svd(WN);
Ww = Dw*Vdw';%WN;%

%==========================================================================
%% Alpha
rng(seed,'twister');
if ~ismember({'Alpha'},who) || isempty(Alpha)
    Alpha = abs(ones(NA,1));%rand(NA,1);%
    Alpha = Alpha/sum(Alpha);
end

%% Init Algorithem Paramters
optionsAlpha = optimset('fmincon');
TypicalAlpha = 0.1*ones(NA,1);
optionsAlpha.Algorithm = 'interior-point';
optionsAlpha.UseParallel = 'always';
optionsAlpha.TolProjCG = 1e-16; 
optionsAlpha.TolProjCGAbs = 1e-20;
optionsAlpha = optimset(optionsAlpha,'LargeScale','on', 'GradObj','on', 'Display','iter', 'TypicalX', TypicalAlpha);
optionsAlpha = optimset(optionsAlpha, 'Hessian',{'lbfgs',30}, 'ScaleProblem', 'obj-and-constr');
optionsAlpha = optimset(optionsAlpha, 'TolX', epsOP, 'TolCon', epsOP, 'TolFun', epsOP,'MaxIter', MaxIter);
%--------------------------------------------------------------------------
%% Main Loop
pause(1);
Start = tic;
funA = @(Alpha)Objfungrad(Alpha, Ww(1:KX,:), epsOP);
A = [];
b = [];
Aeq = ones(1,NA);
beq = 1;
lb = zeros(NA,1);
ub = ones(NA,1);
%----------------------------
[Alpha, J] = fmincon(funA,Alpha,A,b,Aeq,beq,lb,ub,[],optionsAlpha);
toc(Start)
%==========================================================================
%% Show figures
fprintf('Show Figures ... \n');

df = FC(2) - FC(1);
FCE = [FC', FCmax+df:df:Fmax];
lenfc = length(FCE);
indf = find((FCE>=FCmin)&(FCE<=FCmax));
AlphaE = zeros(lenfc,1);

AlphaE(indf) = Alpha/sum(Alpha);
PF = cumsum(AlphaE)/sum(AlphaE);
PF(indf(end):end) = 1.0;

maxa = 0.001*max(AlphaE);
Inda = [100*AlphaE(AlphaE>maxa), FCE(AlphaE>maxa)']
INDA = EvalCluster(Inda)
%--------------------------------------------------------------------------
%%
HFAPA = figure('Name',['Alpha:', FigName]);
[AX,H1,H2] = plotyy(FCE,AlphaE,FCE,PF,'plot');%

xlim(AX(1),[FCmin,Fmax]);
xlim(AX(2),[FCmin,Fmax]);
xlab = [400:50:700];
set(AX(1),'xtick',xlab);
set(AX(2),'xtick',xlab);

ymax = ceil(max(AlphaE)*10)/10;
ymid = ymax/2;
ylab = [0:0.1:max(0.2,ymax)];%[0, ymid, ymax];%
ylim(AX(1),[0,ylab(end)]);
set(AX(1),'ytick',ylab);

ylim(AX(2),[0,1]);
ylab = [0:0.2:1];%
ylim(AX(2),[0,1])
set(AX(2),'ytick',ylab);

%--------------------------------------------------------------------------
set(get(AX(1),'Ylabel'),'String','Population Density {\itp}({\it\theta})')
set(get(AX(2),'Ylabel'),'String','Cumulative Distribution')

set(H1,'LineStyle',':','LineWidth',0.5,'Marker','.','MarkerSize',6,'Color','r')
set(H2,'LineStyle','-','LineWidth',1,'Color','b')%

set(get(AX(1),'Ylabel'),'LineWidth', 1, 'FontSize', FontSize,'Color','r')
set(get(AX(2),'Ylabel'),'LineWidth', 1,'FontSize', FontSize,'Color','b')

box off;
set(AX(1), 'FontSize', FontSize,'XColor', 'k','YColor', 'r');
set(AX(2), 'FontSize', FontSize,'XColor', 'k','YColor', 'b');

xlabel('Wavelength of Maximum Absorbance {\it{\theta}} (nm)');%Peak Sensitivity

% print(HFAPA, '-depsc2','-tiff', '-r600', [FigName, '.eps']);
print(HFAPA, '-dtiff', '-r300', [FigName, '.png']);

%% -----------------------------------------------------------------------
Amax = max(AlphaE);
AlphaP = AlphaE;
AlphaP(AlphaE<0.001*Amax) = 0;
[pks,locs] = findpeaks(AlphaP);%, 'MINPEAKDISTANCE', 50

WI = WI./repmat(max(WI), size(WI,1), 1);
W = W./repmat(max(W), size(W,1), 1);

HFFunFL = figure('Name',['Tuning functions:', figurename]);

COLOR = ['bgrmyc'];%
LenC = length(COLOR);
hold on;
i = 0;
for k = locs'%dk+
    i0 = mod(i,LenC)+1;
    i=i+1;
    color = COLOR(i0);
    plot(F, WI(:,k), ['-',color], 'LineWidth',4)%
end
hold off

ylim([0,1]);

xlim([Fmin,Fmax]);
xlab = [Fmin:50:Fmax];
set(gca,'xtick',xlab);

set(gca,'FontSize', 18,'FontWeight', 'bold','Fontname', 'Times New Roman', 'LineWidth', 1)%

box off;
xlabel('Wavelength (nm)','FontSize', 18,'FontWeight', 'bold','Fontname', 'Times New Roman');
ylabel('Normalized Absorbance','FontSize', 18,'FontWeight', 'bold','Fontname', 'Times New Roman');%

postfix1 = '_t';
% print(HFFunFL, '-depsc2', '-tiff', '-r600', [FigName, postfix1, '.eps']);
print(HFFunFL, '-dtiff', '-r300',[FigName, postfix1, '.png']);

%==========================================================================
%% Save File
filename = ['../results/',figurename,'.mat'];
fprintf('Save File to:  %s \n',filename);
% save(filename, 'Alpha', 'Param', 'optionsAlpha');% 

%==========================================================================
