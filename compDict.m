%init dictionary
Dict = zeros(sparse.dictdimen(1), sparse.dictdimen(2), info.nact+1);
for a = 1:info.nact
    Dict(:,:,a) = normrnd(0,2,sparse.dictdimen);
end
Dict(:,:,info.nact+1) = normrnd(0,2,sparse.dictdimen);

ADim = [(info.nact+1)*sparse.dictdimen(2) 0];
XDim = [0 0];

% init coefficients to 0
for a = 1:info.nact   % update coeff of action a
    idxact = sprintf('a%02d', a);
    load(['Features\',idxact,'.mat']);
    dim = size(feat);
    ADim(2) = ADim(2)+dim(2);
    XDim = [size(feat,1) ADim(2)];
    if exist([info.coefficients, idxact], 'dir')
        rmdir([info.coefficients, idxact], 's');
    end
    mkdir([info.coefficients, idxact]);
    disp(['init action ',idxact, '''s all coefficients block to zero']);
    coeff_common = zeros(sparse.dictdimen(2), dim(2));
    save([info.coefficients, idxact, '\a_s.mat'], 'coeff_common');
    for m = 1:info.nact
        dictindex = sprintf('a%02d', m);
        coeff_specific = zeros(sparse.dictdimen(2), dim(2));
        save([info.coefficients, idxact, '\', dictindex, '.mat'], 'coeff_specific');
    end
end

value = compObjectiveFun(Dict, info, XDim, sparse);

% iteration
for it = 1:sparse.maxiteration
    disp(['iteration times: ', num2str(it), '/', num2str(sparse.maxiteration),...
        ', objective function value=', num2str(value)]);
    
    disp('update the coefficients');
    for a = 1:info.nact   % update coeff of action a
        idxact = sprintf('a%02d', a);
        
        for i = 1:info.nact  % update ai and as
            Da = compDa(Dict, a, dim, info);
            dictindex = sprintf('a%02d', i);
            load([info.coefficients, idxact, '\', dictindex, '.mat']);
            DiAi = Dict(:,:,i) * coeff_specific;  % DiAi
            load([info.coefficients, idxact, '\a_s.mat']);
            DsAs = Dict(:,:,info.nact+1) * coeff_common;  % DsAs
            
            % update ai
            coeff_specific = coeff_specific-sparse.miu*Dict(:,:,i)'*(Da+DiAi+DsAs-2*feat);            
            if i==a
                w = sparse.omiga(1);
            else
                w = sparse.omiga(3);
            end 
            ain2 = norm(coeff_specific,2);
            mlwi = sparse.miu*sparse.lamda*w;
            if  ain2 > mlwi
                coeff_specific = coeff_specific*(ain2-mlwi)/ain2;
            else
                coeff_specific = 0;
            end
            disp(['update action ',idxact, '''s coefficients block a[', num2str(i), ']']);
            save([info.coefficients, idxact, '\', dictindex, '.mat'], 'coeff_specific');
            
            % update as
            coeff_common = coeff_common-sparse.miu*Dict(:,:,info.nact+1)'*(Da+2*DsAs-2*feat);
            asn2 = norm(coeff_common,2);
            mlws = sparse.miu*sparse.lamda*sparse.omiga(2);
            if  asn2 > mlws
                coeff_common = coeff_common*(asn2-mlws)/asn2;
            else
                coeff_common = 0;
            end
            disp(['update action ',idxact, '''s coefficients share block a[s]']);
            save([info.coefficients, idxact, '\a_s.mat'], 'coeff_common');   %%%%%%%%%%%%%%%%???????????????????
        end
    end
    
    disp('update the dictories');
    for Di = 1:info.nact   % update dict Di
        disp(['update the dictories - D',sprintf('%d', Di)]);
        for dj = 1:sparse.dictdimen(2) % update dict dj
            [aj, aij] = getArow(info, Di, dj, info.nact, ADim(2));
            miu1 = 1/norm(aj*aj')^2;
            miu2 = 1/norm(aij*aij')^2;
            
            % comp DA-X
            DA_X = zeros(XDim);
            idxs = 1;
            for i = 1:info.nact
                idxact = sprintf('a%02d', i);
                load(['Features\',idxact,'.mat']);
                dim = size(feat);
                idxe = dim(2)+idxs;
                Dia = compDa(Dict, i, dim, info);
                Dia_x = Dia - feat; 
                DA_X(:,idxs:idxe-1) = Dia_x;
                idxs = idxe;
            end
            tmp1 = miu1*DA_X*aj';
            
            idxact = sprintf('a%02d', Di);
            %comp DiAii
            load([info.coefficients, idxact, '\', idxact, '.mat']);
            DiAii = Dict(:,:,Di)*coeff_specific;
            %comp DsAis
            load([info.coefficients, idxact, '\a_s.mat']);
            DsAis = Dict(:,:,info.nact+1)*coeff_common;

            load(['Features\',idxact,'.mat']);
            tmp2 = miu2*(DiAii+DsAis-feat)*aij';
            
            % comp DDt
            DDt = zeros(size(Dict, 1));
            for didx = 1:info.nact+1
                if didx ~= Di
                    DDt = DDt + Dict(:,:,didx)*Dict(:,:,didx)';
                end
            end
            tmp3 = sparse.yita * DDt;
            D = Dict(:,:,Di);
            djt = D(:,dj);
            djt1 = djt-tmp1-tmp2-tmp3*djt;
            djt1 = djt1/norm(djt1);
            D(:,dj) = djt1;
            Dict(:,:,Di) = D;
        end
    end
    disp(['update the dictories - Ds']);
    for dj = 1:sparse.dictdimen(2)     % update dict Ds
        [aj, ~] = getArow(info, info.nact+1, dj, info.nact, ADim(2));
        miu1 = 1/norm(aj*aj')^2;

        % comp DA-X
        DA_X = zeros(XDim);
        idxs = 1;
        for i = 1:info.nact
            idxact = sprintf('a%02d', i);
            load(['Features\',idxact,'.mat']);
            dim = size(feat);
            idxe = dim(2)+idxs;
            Dia = compDa(Dict, i, dim, info);
            Dia_x = Dia - feat; 
            DA_X(:,idxs:idxe-1) = Dia_x;
            idxs = idxe;
        end
        tmp1 = miu1*DA_X*aj';

        % comp DDt
        DDt = zeros(size(Dict, 1));
        for didx = 1:info.nact
            DDt = DDt + Dict(:,:,didx)*Dict(:,:,didx)';
        end
        tmp2 = sparse.yita * DDt;
        
        tmp3 = zeros(size(Dict,1),1);
        miu2 = 0;
        for i=1:info.nact
            idxact = sprintf('a%02d', i);
            %comp DiAii
            load([info.coefficients, idxact, '\', idxact, '.mat']);
            DiAii = Dict(:,:,i)*coeff_specific;
            %comp DsAis
            load([info.coefficients, idxact, '\a_s.mat']);
            DsAis = Dict(:,:,info.nact+1)*coeff_common;

            load(['Features\',idxact,'.mat']);
            [~, aij] = getArow(info, i, dj, info.nact, ADim(2));
            miu2 = miu2 +norm(aij*aij')^2;
            tmp3 = tmp3 + (DiAii+DsAis-feat)*aij';
        end
        tmp3 = info.nact/miu2*tmp3;
        
        D = Dict(:,:,info.nact+1);
        djt = D(:,dj);
        djt1 = djt-tmp1-tmp2*djt-tmp3;
        djt1 = djt1/norm(djt1);
        D(:,dj) = djt1;
        Dict(:,:,info.nact+1) = D;
    end
    
    value1 = compObjectiveFun(Dict, info, XDim, sparse);
%     if abs(value1-value)/abs(value) <= sparse.threshold
%         disp('has convergenced!');
%         disp(['iterate times: ', num2str(it), ', objective value=', num2str(value1)]);
%         break;
%     else
%         value = value1;
%     end
    value = value1;
end

save('Dict.mat','Dict');

clearvars -except info stip cuboid sparse