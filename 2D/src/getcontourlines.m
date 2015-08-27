function s = getcontourlines(c)

    sz = size(c,2);     % Size of the contour matrix c
    ii = 1;             % Index to keep track of current location
    jj = 1;             % Counter to keep track of % of contour lines

    while ii < sz       % While we haven't exhausted the array
        n = c(2,ii);    % How many points in this contour?
        s(jj).v = c(1,ii);        % Value of the contour
        s(jj).x = c(1,ii+1:ii+n); % X coordinates
        s(jj).y = c(2,ii+1:ii+n); % Y coordinates
        ii = ii + n + 1;          % Skip ahead to next contour line
        jj = jj + 1;              % Increment number of contours
    end

end