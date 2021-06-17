function [A,B,C,D,PI,PI_1] = lfr2pgen(matlfr)
%%
%  File: lfr2pgen.m
%  Directory: 7_ftools/ftools/v12/utilities/lfr_utils
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2020. June 08. (2019b)
%

if ~isa(matlfr,'lfr')
    matlfr = lfr(matlfr);
end

[ny,nu,~,m1_const,m1_row,m1_col] = size(matlfr);
m1 = m1_const + m1_row;
m = nu + m1;

[D,C,B,A,blk] = lfrdata(matlfr);

M = mat2cell(eye(m),m,[nu m1]);

PI = lfr(D,C,M{2},M{1},blk);
PI_1 = lfr(D,C,M{2}(nu+1:end,:),M{1}(nu+1:end,:),blk);

end