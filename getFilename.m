function [filename]=getFilename(path, n1,n2,n3)
    s1 = subr(n1);
    s2 = subr(n2);
    s3 = subr(n3);
    filename = [path, '\a', s1, '_s', s2, '_e', s3, '_sdepth.bin'];
end

function [s]=subr(n)
    if(n<10)
        s = ['0',num2str(n)];
    else
        s = num2str(n);
    end
end