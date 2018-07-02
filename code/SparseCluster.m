function Alpha = SparseCluster(FC, Inda)
    %06/30/2018
    NA = length(FC);
    if ~exist('Inda', 'var')  || isempty(Inda)
        Alpha = abs(ones(NA,1));%rand(NA,1);%
        Alpha = Alpha/sum(Alpha);
        return;
    end
    id = Inda(:, 2);
    val = Inda(:, 1);
    Alpha = zeros(NA,1);
    
    for i = 1:length(id)
        Alpha(FC==id(i)) = val(i);
    end
    Alpha = Alpha/sum(Alpha);
end