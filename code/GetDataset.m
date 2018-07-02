function [X, F] = GetDataset(pdfname, pdfparams)
    % 06/30/2018

    switch pdfname
        case 'CDNS'
            reflectfn = pdfparams.reflectfn;
            Fmin = pdfparams.Fmin;
            Fmax = pdfparams.Fmax;
            load(['../data/', 'ugillums.mat'], 'Lambda', 'illums');
            illums = illums(30, 2:end);
            indx = find((Lambda>=Fmin)&(Lambda<=Fmax));
            F = Lambda(indx);

            if isfield(pdfparams,'r1') 
                a(1) = pdfparams.r1;
                a(2) = pdfparams.r2;
            else
                a(1) = 1;
            end

            Lr = length(reflectfn);
            X = zeros(1, length(Lambda));
            k = 0;
            for n = 1:Lr
                fr = reflectfn{n};
                load(['../data/', fr,'.mat'], 'reflect');
                reflect = reflect(:,2:end);
                spn = size(reflect,1);
                if a(n)<1
                    ind = randsample(spn, round(abs(a(n))*spn));
                    reflect = reflect(ind,:);
                else
                    reflect = repmat(reflect,a(n),1);
                end
                X(k+1:k+size(reflect,1),:) = reflect;
                k = k + size(reflect,1);
            end
            Mr = size(X,1);
            X = (X.*repmat(illums, Mr,1))';
            X = X(indx,:);

        case 'Forest'
            n = pdfparams.typeNo;
            Fmin = pdfparams.Fmin;
            Fmax = pdfparams.Fmax;
            if n == 1 || n == 2 || n == 3
                X = loadForest(n);
            else
                X = [];
                for i = 1:3
                    X = [X, loadForest(i)];
                end
            end
            F = [390:5:850]';
            indx = find((Fmin<=F)&(F<=Fmax));
            X = X(indx,:);
            F = F(indx);

        case 'HINS'
            SampleNum = pdfparams.SampleNum;
            Fmin = pdfparams.Fmin;
            Fmax = pdfparams.Fmax;
            fname = ['../data/','Hypimgs2004.mat'];%
            load(fname,'REFLECT', 'illum_6500', 'paramhys');
            [K, N] = size(REFLECT);
            if SampleNum < N
                rng(pdfparams.seed, 'twister');
                REFLECT = REFLECT(:, randsample(N, SampleNum));
            end
            X = REFLECT.*repmat(illum_6500, 1, SampleNum);
            F = paramhys.Lambda;
            indx = find((Fmin<=F)&(F<=Fmax));
            X = X(indx,:);
            F = F(indx);

        case 'Munsell'
            Fmin = pdfparams.Fmin;
            Fmax = pdfparams.Fmax;
            SampleNum = pdfparams.SampleNum;
            fname = ['../data/munsell380_800_1.mat'];
            load(fname,'munsell');
            [M, N] = size(munsell);
            if N <= SampleNum
                X = munsell;
            else
                rng(pdfparams.seed, 'twister');
                X = munsell(:,randsample(N, SampleNum));
            end
            Lambda = [380:800]';
            indx = find((Fmin<=Lambda)&(Fmax>=Lambda));
            X = X(indx,:);
            F = Lambda(indx);

         case 'Pantone'  
            Fmin = pdfparams.Fmin;
            Fmax = pdfparams.Fmax;
            SampleNum = pdfparams.SampleNum;
            fname = ['../data/pantone.mat'];
            load(fname,'pantone');
            X = pantone;
            F = [380:1:780]';
            N = size(X,2);
            if SampleNum < N
                rng(pdfparams.seed, 'twister');
                X = X(:, randsample(N,SampleNum));
            end
            indx = find((Fmin<=F)&(F<=Fmax));
            X = X(indx,:);
            F = F(indx);
    end
end
%%

function X = loadForest(n)
    if n == 1
        fn = 'pine';
        fname = ['../data/',fn,'.mat'];%
        load(fname);
        X = pine;
    elseif n == 2
        fn = 'spruce';
        fname = ['../data/',fn,'.mat'];%
        load(fname);
        X = spruce;
    elseif n == 3
        fn = 'birch';
        fname = ['../data/',fn,'.mat'];%
        load(fname);
        X = birch;
    else
        X = [];
    end
end
