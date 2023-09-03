function s = movingstd(x,k)
% movingstd: efficient windowed standard deviation of a time series
% usage: s = movingstd(x,k)
%
% Movingstd uses filter to compute the standard deviation, using
% the trick of std = sqrt((sum(x.^2) - n*xbar.^2)/(n-1)).
% Beware that this formula can suffer from numerical problems for
% data which is large in magnitude. Your data is automatically
% centered and scaled to alleviate these problems.
%
% arguments: (input)
%  x   - vector containing time series data
%
%  k   - size of the moving window to use (see windowmode)
%        All windowmodes adjust the window width near the ends of
%        the series as necessary.
%
%        k must be an integer, at least 1 for a 'central' window,
%        and at least 2 for 'forward' or 'backward'
%
%
% arguments: (output)
%  s   - vector containing the windowed standard deviation.
%        length(s) == length(x)

% length of the time series
n = length(x);

% Improve the numerical analysis by subtracting off the series mean
% this has no effect on the standard deviation.
x = x - mean(x);
% scale the data to have unit variance too. will put that
% scale factor back into the result at the end
xstd = std(x);
x = x./xstd;

% we will need the squared elements 
x2 = x.^2;

% sliding std here
A = 1;
B = ones(1,2*k+1);
s = sqrt((filter(B,A,x2) - (filter(B,A,x).^2)*(1/(2*k+1)))/(2*k));
s(k:(n-k)) = s((2*k):end);


% repairs are needed at both ends
s(1:k) = s(k+1);
s(n-k+1:n) = s(n-k);
  
% catch any complex std elements due to numerical precision issues.
% anything that came out with a non-zero imaginary part is
% indistinguishable from zero, so make it so.
s(imag(s) ~= 0) = 0;

% restore the scale factor that was used before to normalize the data
s = s.*xstd;
