function y = pflat(x)
    % Suppose the point is a vertical vector
    y = x./ x(end,:) ; 