for a = 1:info.nact
    for s = 1:info.nsbj
        for e = 1:info.ntms
            load([info.dstippath, getFilename(a, s, e), '_dstip.mat']);
            
            
        end
    end
end

clearvars -except info stip cuboid