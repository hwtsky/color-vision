function INDA = AlphaCluster(Inda)
    id = Inda(:,2);
    L = length(Inda(:,2));
    n = 1;
    j = 1;
    v = Inda(1,1);
    i = Inda(1,2);
    INDA = [];
    for m = 2:L
        if (id(m) - id(m-1) <= 1)
            v = v + Inda(m,1);
            i = i + Inda(m,2);
            n = n + 1;
        else
            INDA(j, 1) = v;
            INDA(j, 2) = i/n;
            j = j + 1;
            n = 1;
            v = Inda(m,1);
            i = Inda(m,2);
        end
        if m == L
            INDA(j, 1) = v;
            INDA(j, 2) = i/n;
        end
    end
end
