function  [aj, aij]  = getArow( info, Di, dj, actioncount, coeffdim)
aj = zeros(1, coeffdim);
idxs = 1;
for i = 1:actioncount
    dictindex = sprintf('a%02d', i);
    if Di ~= info.nact+1
        idxact = sprintf('a%02d', Di);
    else
        idxact = 'a_s';
    end
    load([info.coefficients, dictindex, '\', idxact, '.mat']);
    
    if Di ~= info.nact+1
        idxe = idxs+size(coeff_specific, 2);
        aj(idxs:idxe-1) = coeff_specific(dj,:);
        idxs = idxe;
        if i == Di
            aij = coeff_specific(dj,:);
        end
    else
        idxe = idxs+size(coeff_common, 2);
        aj(idxs:idxe-1) = coeff_common(dj,:);
        idxs = idxe;
        aij = [];
    end
end
end

