function [peaks] = format_rawpeak(raw_peaks)
%FORMAT_RAWPEAK Summary of this function goes here
%   Detailed explanation goes here

peaks = [];
k = 1;
for i=1:length(raw_peaks)-1
    peaks =[peaks;raw_peaks(i,:)];
    peaks =[peaks;raw_peaks(i+1,:)];
end
end

