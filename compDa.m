function Da = compDa(Dict, i, dim, info)
    Da = zeros(dim);
    idxact = sprintf('a%02d', i);
    for a = 1:info.nact   % Da 
        dictindex = sprintf('a%02d', a);
        load([info.coefficients, idxact, '\', dictindex, '.mat']);
        Da = Da + Dict(:,:,a) * coeff_specific;
    end
    load([info.coefficients, idxact, '\a_s.mat']);
    Da = Da + Dict(:,:,info.nact+1) * coeff_common;
end