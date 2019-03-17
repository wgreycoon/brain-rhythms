function [ dataout ] = subtractfirstsamplevalue( datain )

%ensure expected dimensions
if size(datain,2) > size(datain,1)
    datain = datain';
end

%subtract value of first sample from each time series (so all signals start at zero)
d = repmat(datain,1);
for idx = 1:size(datain,2)
    d(:,idx) = datain(:,idx) - datain(1,idx);
end

dataout = d; clear datain d

end