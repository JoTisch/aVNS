%#eml
function[result] = smooth_emlc(data, window)

if length(data) < 2
    result = data;
    return
end

% result = smooth_emlc_worker(data(:), window);

window = 2*ceil(window/2)-1;

result = conv(data, ones(window, 1)/(window));
result = result(ceil(window/2) : end-floor(window/2));

return


function[result] = smooth_emlc_worker(data, window)
window = 2*ceil(window/2)-1;
datalength = length(data);

eins = ones(datalength, 1);

padding = zeros(window-1, 1);
count = [eins; padding];
result = [data; padding];

for i = 1:window-1
    result(i+1:i+datalength) = result(i+1:i+datalength) + data(:,1);
    count(i+1:i+datalength) = count(i+1:i+datalength) + eins;
end

ind = floor(window/2);
count = count(1+ind:end-ind);
result = result(1+ind:end-ind);
result = result ./ count;
return