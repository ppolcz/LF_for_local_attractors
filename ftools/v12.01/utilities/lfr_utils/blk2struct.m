function [ret] = blk2struct(blk)
%%
%  File: blk2struct.m
%  Directory: 7_ftools/ftools/v12/utilities/lfr_utils
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2020. June 08. (2019b)
%

if isa(blk,'plfr') || isa(blk,'lfr') 
    blk = blk.blk;
end

names = blk.names;
desc = blk.desc;

if isempty(names)
    ret = struct;
    return
end

if strcmp(names{1},'1')
    names{1} = 'CONST__';
end

if size(desc,1) < 13
    desc(12:13,:) = [-Inf;Inf]*ones(1,size(desc,2));
end

ret = cell2struct(num2cell(desc,1),names,2);

end


function test1_lfr_with_constant_block
%%
lfrs a b c d

G = 1/a + b*c*d + 1/(a+b+d+10);

s = blk2struct(G)

end