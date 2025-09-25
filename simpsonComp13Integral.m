function cost = simpsonComp13Integral(t, y)
    tic 
    N = length(t) - 1;
    if mod(N, 2) ~= 0
        error('Simpson''s rule requires an even number of intervals');
    end
    h = (t(end) - t(1)) / N;
    cost = y(1) + y(end) + 4*sum(y(2:2:end-1)) + 2*sum(y(3:2:end-2));
    cost = cost * (h/3);
    toc
end