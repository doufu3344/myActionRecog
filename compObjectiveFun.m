function value = compObjectiveFun( Dict, info, XDim, sparse)
value = 0;
valuetmp = 0;

DA_X = zeros(XDim);
idxs = 1;
for i = 1:info.nact
    idxact = sprintf('a%02d', i);
    load(['Features\',idxact,'.mat']);
    
    load([info.coefficients, idxact, '\', idxact, '.mat']);
    DiAii = Dict(:,:,i)*coeff_specific;
    %comp DsAis
    load([info.coefficients, idxact, '\a_s.mat']);
    DsAis = Dict(:,:,info.nact+1)*coeff_common;
    value = value+norm(feat-DiAii-DsAis, 'fro')^2;
    
    D = zeros(size(Dict, 2));
    for tmp = i+1:info.nact
        D = D + Dict(:,:,i)'*Dict(:,:,tmp);
        value = value+sparse.yita*norm(D, 'fro')^2;
    end

    dim = size(feat);
    idxe = dim(2)+idxs;
    Dia = compDa(Dict, i, dim, info);
    Dia_x = feat - Dia; 
    DA_X(:,idxs:idxe-1) = Dia_x;
    idxs = idxe;
    
    for cof = 1:info.nact
        coffidx = sprintf('a%02d', cof);
        load([info.coefficients, idxact, '\', coffidx, '.mat']);
        if cof ==i
            valuetmp = valuetmp+ norm(coeff_specific, 2)*sparse.omiga(1);
        else
            valuetmp = valuetmp+ norm(coeff_specific, 2)*sparse.omiga(3);
        end
    end
    load([info.coefficients, idxact, '\a_s.mat']);
    valuetmp = valuetmp+ norm(coeff_common, 2)*sparse.omiga(2);
end
value = value+norm(DA_X, 'fro')^2;
value = value+valuetmp*sparse.lamda;
end

