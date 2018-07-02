function [J, G] = Objfungrad(Alpha, W, EPS)
    % 06/30/2018 

    [KX, NA] = size(W);
    if ~exist('EPS', 'var')
        EPS = eps;
    end
    J = 0;
    GA = 0;
    Wg = W.*repmat(sqrt(Alpha)', KX,1);
    L = Wg*Wg' + diag(ones(KX,1));
    [LG, pl] = chol(L,'lower');
    J1 = 2*sum(log(diag(LG)+EPS));
    J = J + J1;
    %------------------------------
    flag  = (nargout > 1);
    if flag
        ga1 = sum(L\W.*W)';
        GA = GA + ga1;
    end
    J = -0.5*J;
    if flag
        G = (-0.5)*GA;
    end
end

