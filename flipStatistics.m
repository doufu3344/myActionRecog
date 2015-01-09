% calclcuate the number of flip of the signal.
% output: nflip --number of flip
%         avgdura --average duration of the flip
function [nflip, avgdura] = flipStatistics( signal )
    signal = signal-mean(signal);
    signal(abs(signal)<50) = 0;
    s = sign(signal);

    t = filter([1 1], 1, s);
    nflip = (length(s)-length(find(t)));

    positive = sum(s>0);
    negative = sum(s<1);

    if positive > negative
        s=-s;
    end
    duration = sum(s>-1);
    avgdura=duration/(nflip+1e-9);
end

