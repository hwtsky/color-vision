function [J, G] = Objfungrad(Alpha, W, EPS)
    % 06/30/2018 
    KX = size(W, 1);
    if ~exist('EPS', 'var')
        EPS = eps;
    end
    Wg = W.*repmat(sqrt(Alpha)', KX,1);
    L = Wg*Wg' + diag(ones(KX,1));
    [LG, pl] = chol(L,'lower');
    J = -sum(log(diag(LG)+EPS));
    if nargout > 1
        G = -0.5*sum(L\W.*W)';
    end
end