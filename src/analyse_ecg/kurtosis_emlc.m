%#eml
function k = kurtosis_emlc(x,flag,dim)
%KURTOSIS Kurtosis.
%   K = KURTOSIS(X) returns the sample kurtosis of the values in X.  For a
%   vector input, K is the fourth central moment of X, divided by fourth
%   power of its standard deviation.  For a matrix input, K is a row vector
%   containing the sample kurtosis of each column of X.  For N-D arrays,
%   KURTOSIS operates along the first non-singleton dimension.
%
%   KURTOSIS(X,0) adjusts the kurtosis for bias.  KURTOSIS(X,1) is the same
%   as KURTOSIS(X), and does not adjust for bias.
%
%   KURTOSIS(X,FLAG,DIM) takes the kurtosis along dimension DIM of X.
%
%   KURTOSIS treats NaNs as missing values, and removes them.
%
%   See also MEAN, MOMENT, STD, VAR, SKEWNESS.

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision: 1.1.8.1 $  $Date: 2010/03/16 00:15:06 $

if nargin < 2 || isempty(flag)
    flag = 1;
end
if nargin < 3 || isempty(dim)
    % The output size for [] is a special case, handle it here.
    if isequal(x,[]), k = NaN; return; end;

    % Figure out which dimension nanmean will work along.
    dim = find(size(x) ~= 1, 1);
    if isempty(dim), dim = 1; end
end


% Center X, compute its fourth and second moments, and compute the
% uncorrected kurtosis.

% mittelwert aus vektorsubtraktion herausnehmen -> MASSIVES
% speedup (weil mean sonst für jedes element im vektor neu berechnet wird)
n = size(x,dim);
mittelwert = sum(x,dim)/n;
x = x - mittelwert;
x2 = x.^2;
s2 = sum(x2,dim)/n; % this is the biased variance estimator
m4 = sum(x2.^2,dim)/n;
k = m4 ./ s2.^2;

% Bias correct the kurtosis.
% removed for speedup
% if flag == 0
%     n = sum(~isnan(x),dim);
%     n(n<4) = NaN; % bias correction is not defined for n < 4.
%     k = ((n+1).*k - 3.*(n-1)) .* (n-1)./((n-2).*(n-3)) + 3;
% end





function m = nanmean(x,dim)
%NANMEAN Mean value, ignoring NaNs.
%   M = NANMEAN(X) returns the sample mean of X, treating NaNs as missing
%   values.  For vector input, M is the mean value of the non-NaN elements
%   in X.  For matrix input, M is a row vector containing the mean value of
%   non-NaN elements in each column.  For N-D arrays, NANMEAN operates
%   along the first non-singleton dimension.
%
%   NANMEAN(X,DIM) takes the mean along dimension DIM of X.
%
%   See also MEAN, NANMEDIAN, NANSTD, NANVAR, NANMIN, NANMAX, NANSUM.

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision: 1.1.8.1 $  $Date: 2010/03/16 00:15:50 $


%%%% Replaced by standard mean without sanity check for speedup
m = sum(x,dim)/size(x,dim);
return

% Find NaNs and set them to zero
nans = isnan(x);
x(nans) = 0;

if nargin == 1 % let sum deal with figuring out which dimension to use
    % Count up non-NaNs.
    n = sum(~nans);
    n(n==0) = NaN; % prevent divideByZero warnings
    % Sum up non-NaNs, and divide by the number of non-NaNs.
    m = sum(x) ./ n;
else
    % Count up non-NaNs.
    n = sum(~nans,dim);
    n(n==0) = NaN; % prevent divideByZero warnings
    % Sum up non-NaNs, and divide by the number of non-NaNs.
    m = sum(x,dim) ./ n;
end

