for a = 1:info.nact
    for s = 1:info.nsbj
        for e = 1:info.ntms
            filename = [info.vidpath, getFilename(a, s, e), '_sdepth.bin'];
            disp(['computing normals of video: ', getFilename(a, s, e), '_sdepth.bin']);
            depth = readDepthBin(filename);

            % some missed videos
            if isempty(depth)
                continue;
            end
            
            % compute derivatives of depth sequence
            [nrows, ncols, nfrms] = size(depth);
            dx = zeros(nrows, ncols, nfrms - 1);
            dy = zeros(nrows, ncols, nfrms - 1);
            dt = zeros(nrows, ncols, nfrms - 1);
            mag = zeros(nrows, ncols, nfrms - 1);
            
            for f = 1:nfrms-1
                % smooth
                frame1 = medfilt2(depth(:, :, f), [5, 5]);
                frame2 = medfilt2(depth(:, :, f + 1), [5, 5]);
                
                % derivatives along x/y/t
                [dx(:, :, f), dy(:, :, f)] = gradient(frame1);
                dt(:, :, f) = frame2 - frame1;
            
                % normalize
                reg = sqrt(dx(:, :, f).^2 + dy(:, :, f).^2 + dt(:, :, f).^2);
                dx(:, :, f) = dx(:, :, f) ./ reg;
                dy(:, :, f) = dy(:, :, f) ./ reg;
                dt(:, :, f) = dt(:, :, f) ./ reg;
                mag(:, :, f) = reg;
                
                dx(isinf(dx)) = 0; dx(isnan(dx)) = 0;
                dy(isinf(dy)) = 0; dy(isnan(dy)) = 0;
                dt(isinf(dt)) = 0; dt(isnan(dt)) = 0;
                mag(isinf(mag)) = 0; mag(isnan(mag)) = 0;
            end
            
            % save normals and masks
            normalName = [info.normalpath, getFilename(a, s, e), '_norm.mat'];
            save(normalName, 'dx', 'dy', 'dt', 'mag');
        end
    end
end

clearvars -except info stip cuboid