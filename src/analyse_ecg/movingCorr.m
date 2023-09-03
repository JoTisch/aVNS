function Cor = movingCorr(template, y)
k = length(template);
n = length(y);

if (n<k)
    Cor = NaN(n,1);
else
    A = 1;
    B = ones(1,k);
    m = filter(B,A,template)/k;
    x = template - m(end);
    x2 = x.^2;
    y2 = y.^2;
    sxy = filter(x(end:-1:1), 1, y);
    Stdx = sqrt((filter(B,A,x2) - (filter(B,A,x).^2)*(1/k))/(k-1));
    Stdy = sqrt((filter(B,A,y2) - (filter(B,A,y).^2)*(1/k))/(k-1));
    Cor = sxy ./ ((k-1)*Stdx(end)*Stdy);
    Cor = Cor(k:end);
end
