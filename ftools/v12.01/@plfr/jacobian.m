function JLFR = jacobian(pLFR,p_cell,subsvars)
%%

JLFR = zeros(pLFR.ny,0);

for i = 1:numel(p_cell)
    try
        JLFR = [JLFR , diff(pLFR.lfrtbx_obj,p_cell{i}.blk.names{1}) ];
    catch e
        if strfind(e.message, 'does not exist in lfr-object')
            continue
        else
            getReport(e)
            break;
        end
    end
end

JLFR = plfr(JLFR);

% dp_plfr = plfr(vertcat(dp_cell{:}));
% 
% if nargin < 4
%     subsvars = plfr.var_helper__(pLFR,dp_plfr);
% end
% 
% dLFR = dLFR.set_vars(subsvars(:));

end