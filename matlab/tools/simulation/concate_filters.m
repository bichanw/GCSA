function B = concate_filters(flts)
% concatenate filters for eigenvalue calculation to prevent blowing up
% start concatenating from the second filter because I wrote this function thinking that the first filter is for input
B = cat(2,flts(2:end).wfilt);
ny = size(B,1);
for ii = 1:(numel(flts)-2)
    tmp = [];
    for jj = 1:(numel(flts)-1)
        if jj == ii
            tmp = cat(2,tmp,eye(ny));
        else
            tmp = cat(2,tmp,zeros(ny));
        end
    end
    B = cat(1,B,tmp);  
end

end