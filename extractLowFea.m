for a = 1:info.nact
    for s = 1:info.nsbj
        for e = 1:info.ntms
            disp(['computing feature of video: ', getFilename(a, s, e), '_sdepth.bin']);
            load([info.dstippath, getFilename(a, s, e), '_dstip.mat']);           
            npoints = length(dstip);
            
            load([info.normalpath, getFilename(a, s, e), '_norm.mat']);

            depth = readDepthBin([info.vidpath, getFilename(a, s, e), '_sdepth.bin']);
            [nrows, ncols, nfrm] = size(depth);
            
            feat = zeros(npoints, 4*cuboid.xpoints*cuboid.ypoints*cuboid.temporalscale);
            for i = 1:npoints
                xi = dstip(i, 1); yi = dstip(i, 2); ti = dstip(i, 3);
                di = depth(xi,yi,ti);
                if ti>floor(cuboid.temporalscale/2)
                    di = di + sum(depth(xi,yi,ti-floor(cuboid.temporalscale/2):ti-1));
                end
                if ti<=(nfrm-floor(cuboid.temporalscale/2))
                    di = di + sum(depth(xi,yi,ti+1:ti+floor(cuboid.temporalscale/2)));
                end
                di = di/cuboid.temporalscale;
                if di ~= 0 && di < cuboid.spatioscale
                    step = floor(cuboid.spatioscale/di);
                else
                    step = 1;
                end
                
                [nx, ny, nt] = size(dx);
                xedge = floor(cuboid.xpoints/2)*step;
                yedge = floor(cuboid.ypoints/2)*step;
                tedge = floor(cuboid.temporalscale/2);
                DX = zeros(nx+xedge*2, ny+yedge*2, nt+tedge*2);
                DX(xedge+1:xedge+nx, yedge+1:yedge+ny, tedge+1:tedge+nt) = dx;
                DY = zeros(nx+xedge*2, ny+yedge*2, nt+tedge*2);
                DY(xedge+1:xedge+nx, yedge+1:yedge+ny, tedge+1:tedge+nt) = dy;
                DT = zeros(nx+xedge*2, ny+yedge*2, nt+tedge*2);
                DT(xedge+1:xedge+nx, yedge+1:yedge+ny, tedge+1:tedge+nt) = dt;
                MAG = zeros(nx+xedge*2, ny+yedge*2, nt+tedge*2);
                MAG(xedge+1:xedge+nx, yedge+1:yedge+ny, tedge+1:tedge+nt) = mag;
                
                xs = xi; xe = xi+floor(cuboid.xpoints/2)*step*2;
                ys = yi; ye = yi+floor(cuboid.ypoints/2)*step*2;
                ts = ti; te = ti+floor(cuboid.temporalscale/2)*2;
                
                featix = DX( xs:step:xe, ys:step:ye,ts:te );
                featiy = DY( xs:step:xe, ys:step:ye,ts:te );
                featit = DT( xs:step:xe, ys:step:ye,ts:te );
                featim = MAG( xs:step:xe, ys:step:ye,ts:te );
                clear DX DY DT MAG;
                feati = [featix(:); featiy(:); featit(:); featim(:)];
                feat(i,:) = feati;
            end
            featName = [info.featurespath, getFilename(a, s, e), '_feat.mat'];
            save(featName, 'feat');           
        end
    end
end

clearvars -except info stip cuboid sparse