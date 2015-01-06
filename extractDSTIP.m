for a = 6:info.nact
    for s = 1:info.nsbj
        for e = 1:info.ntms
            % Read depth sequence
            filename = getFilename(info.vidpath, a, s, e);  
            depth_ori = readDepthBin(filename);
            [nrows, ncols, nfrm] = size(depth_ori);
            depth_filtered = zeros(nrows, ncols, nfrm);
            if nfrm < 2*stip.gabor_scale+1
                error('Gabor filter scale is too large!');
            end
            
            % Apply gaussian filter
            h = fspecial('gaussian', [stip.gauss_size, stip.gauss_size], stip.gauss_variance);
            for f = 1:size(depth_ori, 3)
                depth_filtered(:,:,f) = conv2( depth_ori(:,:,f), h, 'same');
            end
            
            % Apply gabor filter
            tao = stip.gabor_scale;
            hev = zeros(1,2*tao+1);
            hod = zeros(1,2*tao+1);
            w = stip.gabor_Ormig;
            t = -tao:tao;
            hev(tao+t+1) = -cos(2*pi*t*w).*exp(-t.^2/tao^2);
            hod(tao+t+1) = -sin(2*pi*t*w).*exp(-t.^2/tao^2);
            hev = hev-mean(hev);
            hod = hod-mean(hev);
            Rev = convn(depth_filtered, shiftdim(hev,-1), 'same');
            Rod = convn(depth_filtered, shiftdim(hod,-1), 'same');
            depth_filtered = Rev+Rod;
            clear Rev Rod;
            
            % Noise suppressor
            if stip.Noisesup_all == 0 % use time interval
                frame_flip = zeros(nrows, ncols, nfrm);
                for i = 1:nrows
                    for j = 1:ncols
                        for t = tao+1:nfrm-tao-1
                            signal = depth_ori(i, j, t-tao:t+tao);
                            signal = signal(:);
                            signal = signal-mean(signal);
                            [n,ad]= flipStatistics(signal);
                            frame_flip(i,j,t)=ad;
                        end
                    end
                end
                depth_filtered=depth_filtered.*frame_flip;
            else
                frame_flip=zeros(nrows,ncols);
                for i = 1:nrows
                    for j = 1:ncols
                        signal = depth_ori(i,j,:);
                        signal = signal(:);
                        [n,ad]= flipStatistics(signal);
                        frame_flip(i,j)=ad;
                    end
                end
                frame_flip(frame_flip < stip.noise_thresh) = 0;
                depth_filtered=depth_filtered.*repmat(frame_flip,[1 1 nfrm]);
            end           
            % extract maxmum
            depth_filtered=abs(smooth3(depth_filtered(:,:,:),'gaussian',[7 7 3],1)); % smooth the data
            depth_filtered(:,:,1:tao)=0;
            depth_filtered(:,:,end-tao:end)=0;
            M = imregionalmax(depth_filtered);  % find local maximum
            S = sort(depth_filtered(M));
            if length(S) > stip.npoints  % take the largest n points
                thresh=S(end - stip.npoints+1);
                M(depth_filtered < thresh)=0;
            end
            
            LinearIndex=find(M==1);   
            [I1,I2,I3]=ind2sub(size(M),LinearIndex);
            %plot(I2,I1,'r*');  %plot(I2,nrows-I1,'*');
            %axis([0 ncols 0 nrows]);
            %clear LinearIndex M;
            
            for tt = 1:nfrm
                h = figure;
                imagesc(depth_ori(:,:,tt));
                hold on;
                plot(I2,I1,'y*');  %plot(I2,nrows-I1,'*');
                axis([0 ncols 0 nrows]);
                saveas(h,['pic\f_', num2str(tt),'.jpg'],'jpg');
                close(h);
            end
        end
    end
end