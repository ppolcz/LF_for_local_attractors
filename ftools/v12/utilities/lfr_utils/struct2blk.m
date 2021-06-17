function [ret] = struct2blk(struct)
%%
%  File: struct2blk.m
%  Directory: 7_ftools/ftools/v12/utilities/lfr_utils
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2020. June 08. (2019b)
%

if isempty(struct)
    ret.names = cell(1,0);
    ret.desc = zeros(13,0);
end

names = fieldnames(struct);

if strcmp(names{1},'CONST__')
    names{1} = '1';
end


desc = struct2cell(struct);
desc = [ desc{:} ];


ret.desc = desc;
ret.names = names(:).';

end


function test1_lfr_with_constant_block
%% syslfr = lfrbounds(syslfr) [LENYEGEBEN]

G_reset(1)

lfrs a b c d

G = 1/a + b*c*d + 1/(a+b+d+10);

s = blk2struct(G);

blk = struct2blk(s);

[A,B,C,D,blk_eredeti] = lfrdata(G);

G_uj = lfr(A,B,C,D,blk);

display(blk_eredeti)
display(blk)

pcz_symzero(lfr2sym(G) - lfr2sym(G_uj),'The two LFRs are practically the same')


end