function saveStipPic( depth, I1, I2)
h = figure;
for f = 1:size(depth, 3)
    imagesc(depth(:,:,f));
    hold on;
    plot(I2,I1,'y*');  %plot(I2,nrows-I1,'*');
    axis([0 size(depth,2) 0 size(depth,1)]);
    saveas(h, ['pic\f_', num2str(f),'.jpg'], 'jpg');    
    clf(h);
end
end

