function [ret] = lfrbounds(syslfr,w)
%%
%  File: lfrbounds.m
%  Directory: 7_ftools/ftools/v12/utilities/lfr_utils
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2020. June 08. (2019b)
%

if isa(syslfr,'plfr')
    syslfr = syslfr.lfrtbx_obj;
end
blk = syslfr.blk;

if isa(w,'plfr')
    w = w.lfrtbx_obj;
end

if nargin == 1
    blk = struct2blk(blk2struct(syslfr));
elseif nargin == 2
    
    % In blk_tmp, the block dimension arguments do not correspond to the
    % dimensions of syslfr. However, the bound information are updated.
    blk_tmp = struct2blk(...
        parsepropval('ignore',blk2struct(syslfr),blk2struct(w)));
    
    % Update the bound information in the original blk structure.
    blk.desc(9:13,:) = blk_tmp.desc(9:13,:);
end

[A,B,C,D] = lfrdata(syslfr);
ret = lfr(A,B,C,D,blk);

end