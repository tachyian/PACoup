function [ lfp_mean ] = get_lfp_mean( lfp_array,  m1array_to_remove)
%GET_LFP_MEAN Get mean LFP from a 10 by 10 array of LFPs
%   Detailed explanation goes here

lfp_mean = zeros( size(lfp_array{1}{1} ) );
if ~isempty(lfp_array)
    n_mean = 0;
    for row=1:10
        for col=1:10
            if ~( ( row == 1 && col == 1 ) || ( row == 1 && col == 10 ) ...
                    || ( row == 10 && col == 1 ) || ( row == 10 && col == 10 ) ) ...
                    && (isempty (strmatch([row col], m1array_to_remove)))
                lfp_mean = lfp_mean + lfp_array{row}{col};
                n_mean = n_mean +1;
            end
        end
    end
    lfp_mean = lfp_mean / n_mean;
end


end

