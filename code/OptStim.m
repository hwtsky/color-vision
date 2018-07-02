%% OptStim 06/30/2018
%==========================================================================
clear all
close all
clc

%% Init Parameters
fprintf('\nStart OptStim ...\n');
%--------------------------------------------------------------------------
pdfname = 'CDNS'; % 'CDNS' 'Forest' 'HINS' 'Munsell'
typeNo = 0; %0;% 1;% 2;% 3;%
%------------------------------------- 
fileNo = '3-1';% '3-2';%
postfix = 'a';
%-------------------------------------
DF = 1;
DFC = 1;%0.5;%
Fmin = 380;% >= 380
Fmax = 700;
FCmin = 350;%Fmin;%
FCmax = 565;% 565 570 580 600 700
F = [Fmin:DF:Fmax]';
FC = [FCmin:DFC:FCmax]';
%% Alpha
KK = 5;
if ~ismember({'Alpha'},who) || isempty(Alpha)
    Inda = [];%[[55, 45]', [429, 555]'];%[[5, 30, 65]', [420, 535, 565]'];%[[25, 25, 25, 25]', [372, 449, 502, 563]'];%[];%[];%[[1, 6, 12]', [430, 540, 565]'];%
    Alpha = SparseCluster(FC, Inda);
end

%% ------------------------------------------------------------------------
figurename = ['Fig', num2str(fileNo), '_', postfix]; 
FigName = ['../results/figures/', figurename];

FontSize = 10;
LineWidth = 0.5;
seed = typeNo + 1;
SampleNum = 10000;
paramGA.Dpeak = 0.3;
paramGA.Trans = 1;
paramGA.FgN = 0;
paramGA.r = 0;
paramGA.Lc = 1.0;
paramGA.Mc = 1.0;
paramGA.isPoly = 0;
NA = length(FC);
[W, paramGA, WI] = Wpig(FC, F, paramGA);

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
    
    %% ----------------------------------------------------------------------
%     [NX, SampleNum] = size(X);
%     LS = 1:SampleNum;
%     [x2, idx] = sort(sum(X.*X));
%     X = X(:,idx);

%     HX = figure('Name',['X'],'PaperPosition', [1.5, 7.5, 5.0, 3]);
%     axes('Position',[0.15 0.15 0.7 0.7]);
%     imagesc(F, LS, X');
% %     surf(LS, F,(X(1:NX,:)+eps),'edgecolor','none'); view(0,90); axis tight;%Flg2
%     box on;
%     set(gca,'XDir','normal')
%     set(gca,'YDir','normal')
%     set(gca, 'color', 'blue');
%     colormap(mycolor(256));
%     bartic = [min(X(:)), max(X(:))]';%%[0:0.2:1]';%
%     barlabel = {'min', 'max'};%bartic;%
%     hcb = colorbar('Position',[0.87 0.35 0.03 0.35],'ytick',bartic,'YTickLabel',barlabel);%,'TickDir','in'
%     set(hcb,'YTickMode','manual','FontSize', FontSize-2, 'Fontname', 'Arial','LineWidth', LineWidth,'color','k'); %'FontWeight', 'bold',
%     
%     ylabel('Reflectance spectrum samples (#)','FontSize',FontSize,'Fontname', 'Arial','color','k'); 
%     xlabel('Wavelength (nm)','FontSize',FontSize,'Fontname', 'Arial','color','k');
%     xtick = [Fmin:50:Fmax];
%     ytick = [0:50:SampleNum];
%     set(gca,'xtick',xtick,'XTickLabel',xtick, 'TickDir', 'out','FontSize', FontSize-2, 'Fontname', 'Arial','LineWidth', LineWidth);% 
%     set(gca,'TickDir', 'out','FontSize', FontSize-2, 'Fontname', 'Arial','LineWidth', LineWidth);% 'ytick',ytick,'YTickLabel',ytick,'TickDir', 

%     print(HX, '-depsc2 ', '-r600',[FigName, '_', 'X','.eps']);
    
    %% --------------------------------------------------------------------
%     HFXm = figure('Name',['Mean of Inputs']);
%     Xmean = mean(X,2);
%     Xmax = max(Xmean);
%     hplot = plot(F, Xmean,'r','LineWidth',3);%/Xmax
%     haxis = get( hplot, 'Parent' ); 
% 
%     xlim(gca,[Fmin,Fmax]);
%     xtick = [Fmin:50:Fmax];
%     set(gca,'xtick',xtick);
% 
%     ymax = ceil(Xmax*1000)/1000;%1;%
%     ymid = ymax/2;
%     ytick0 = [0, ymid, ymax];%
%     ylim([0,ytick0(end)]);
%     set(gca,'ytick',ytick0);
%     
%     set(gca, 'FontSize', FontSize)%, 'FontWeight', 'bold','Fontname', 'Times New Roman'
%     
%     box off;
%     xlabel('Wavelength (nm)');
%     ylabel('Mean of Stimulus Spectra');%Amplitude %of Stimulus Spectra
 
%     print(HFXm, '-depsc2 ', '-r600',[FigName, '_', 'Xmean','.eps']);
    %-----------------------------
end
%% SVD
[NX, SampleNum] = size(X);
xv = mean(X, 2);
X = X - repmat(xv, 1, SampleNum);
[U, D, Vd] = svd(X);%
Vd = [];
% -----------------------------------------------------
B = W*diag(sqrt(Alpha));
[Ub, Db, Vb] = svd(B);%

%%
Uf = Ub;%  U; %

if sum(Uf(:,1)) > 0
    Uf = -Uf;
end

HFPCS = figure('Name',['Uf:', num2str(KK)]);
hold on

xlim([Fmin,Fmax]);
plot(F,Uf(:,1),'r-','LineWidth',2)
plot(F,Uf(:,2),'m-','LineWidth',2)

if KK <= 2
    HLEG = legend('1','2','Location','SouthWest');
    set(HLEG, 'FontSize', FontSize);
end
if KK >= 3
    plot(F,Uf(:,3),'b-','LineWidth',2)
    HLEG = legend('1','2','3','Location','SouthWest');
    set(HLEG, 'FontSize', FontSize);
end

if KK>=4
    semilogx(F,Uf(:,4),'g-','LineWidth',2)
    HLEG=legend('1','2','3','4','Location','SouthWest');
    set(HLEG, 'FontSize', FontSize);
end

if KK >= 5
    semilogx(F,Uf(:,5),'c-','LineWidth',2)
%     ylim([min(U(:,5)),max(U(:,5))]);
    HLEG=legend('1','2','3','4','5','Location','NorthEast');
    set(HLEG, 'FontSize', FontSize);
end

Uk = Uf(:,1:KK);
Ymin = -0.2;%min(Uk(:));
Ymax = 0.2;%max(Uk(:));
ymax = ceil(Ymax*10)/10;
ymid = ymax/2;
ymin = floor(Ymin*10)/10;
ylab = [-max(abs([ymin,ymax])),0, max(abs([ymin,ymax]))];%
set(gca,'ytick',ylab);
ylim([ylab(1),ylab(end)]);

xlim(gca,[Fmin,Fmax]);
xlab = [Fmin:50:Fmax];
set(gca,'xtick',xlab);

set(gca, 'FontSize', FontSize)

box off;
xlabel('Wavelength (nm)');
ylabel('Principal Components');

% print(HFPCS, '-depsc2 ', '-r600',[FigName, '.eps']);
print(HFPCS, '-dtiff', '-r300', [FigName, '.png']);

%% Show figures
fprintf('Show Figures ... \n');

df = FC(2) - FC(1);
FCE = [FC', FCmax+df:df:Fmax];
lenfc = length(FCE);
indf = find((FCE>=FCmin)&(FCE<=FCmax));
AlphaE = zeros(lenfc,1);%

AlphaE(indf) = Alpha/sum(Alpha);
PF = cumsum(AlphaE)/sum(AlphaE);
PF(indf(end):end) = 1.0;
%--------------------------------------------------------------------------
%%
% HFAPA = figure('Name',['Alpha:',fn]);
% [AX,H1,H2] = plotyy(FCE,AlphaE,FCE,PF,'plot');%'semilogy',
% 
% xlim(AX(1),[FCmin,Fmax]);
% xlim(AX(2),[FCmin,Fmax]);
% xlab = [400:50:700];
% set(AX(1),'xtick',xlab);
% set(AX(2),'xtick',xlab);
% 
% ymax = ceil(max(AlphaE)*10)/10;
% ymid = ymax/2;
% ylab = [0:0.1:max(0.2,ymax)];%[0, ymid, ymax];%
% ylim(AX(1),[0,ylab(end)]);
% set(AX(1),'ytick',ylab);
% 
% ylim(AX(2),[0,1]);
% ylab = [0:0.2:1];%
% ylim(AX(2),[0,1])
% set(AX(2),'ytick',ylab);
% 
% %--------------------------------------------------------------------------
% set(get(AX(1),'Ylabel'),'String','Population Density {\itp}({\it\theta})')
% set(get(AX(2),'Ylabel'),'String','Cumulative Distribution')
% 
% set(H1,'LineStyle',':','LineWidth',0.5,'Marker','.','MarkerSize',6,'Color','r')%
% set(H2,'LineStyle','-','LineWidth',1,'Color','b')%
% 
% set(get(AX(1),'Ylabel'),'LineWidth', 1, 'FontSize', FontSize,'Color','r')%'Fontname', 'Times New Roman','FontWeight', 'bold', 
% set(get(AX(2),'Ylabel'),'LineWidth', 1,'FontSize', FontSize,'Color','b')%'Fontname', 'Times New Roman','FontWeight', 'bold', 
% 
% box off;
% set(AX(1), 'FontSize', FontSize,'XColor', 'k','YColor', 'r');%'Fontname', 'Times New Roman', 'FontWeight', 'bold','LineWidth', 1,
% set(AX(2), 'FontSize', FontSize,'XColor', 'k','YColor', 'b');%'Fontname', 'Times New Roman', 'FontWeight', 'bold','LineWidth', 1,
% 
% xlabel('Wavelength of Maximum Absorbance {\it{\theta}} (nm)');%Peak Sensitivity

% print(HFAPA, '-depsc2','-tiff', '-r600', [FigName, '.eps']);
% print(HFAPA, '-dtiff','-r300',[FigName, '.png']);

%% -----------------------------------------------------------------------
% Amax = max(AlphaE);
% AlphaP = AlphaE;
% AlphaP(AlphaE<0.001*Amax) = 0;
% [pks,locs] = findpeaks(AlphaP);%, 'MINPEAKDISTANCE', 50
% 
% WI = WI./repmat(max(WI), size(WI,1), 1);
% W = W./repmat(max(W), size(W,1), 1);
% 
% HFFunFL = figure('Name',['Tuning functions:', figurename]);
% 
% COLOR = ['bgrmyc'];%
% LenC = length(COLOR);
% hold on;
% i = 0;
% for k = locs'%dk+
%     i0 = mod(i,LenC)+1;
%     i=i+1;
%     color = COLOR(i0);
%     plot(F, WI(:,k), ['-',color], 'LineWidth',4)%
% end
% hold off
% 
% ylim([0,1]);
% 
% xlim([Fmin,Fmax]);
% xlab = [Fmin:50:Fmax];
% set(gca,'xtick', xlab);
% 
% set(gca,'FontSize', 18,'FontWeight', 'bold','Fontname', 'Times New Roman', 'LineWidth', 1)%
% 
% box off;
% xlabel('Wavelength (nm)','FontSize', 18,'FontWeight', 'bold','Fontname', 'Times New Roman');
% ylabel('Normalized Absorbance','FontSize', 18,'FontWeight', 'bold','Fontname', 'Times New Roman');%
% postfix1 = '_t';
% print(HFFunFL, '-depsc2','-tiff', '-r600', [FigName, postfix1, '.eps']);
% print(HFFunFL, '-dtiff','-r300',[FigName, postfix1, '.png']);

%==========================================================================