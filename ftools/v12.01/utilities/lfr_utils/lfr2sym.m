function [pLFR_sym, PI_sym] = lfr2sym(matlfr)
%%
%  File: lfr2sym.m
%  Directory: 7_ftools/ftools/v12/utilities/lfr_utils
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2020. June 08. (2019b)
%

if isa(matlfr,'plfr')
    [pLFR_sym, PI_sym] = sym(matlfr);
    return;
end

[A,B,C,D,Delta] = data(matlfr);  

[ny,nu,nx,m1_const,m1_row,m1_col] = size(matlfr);

I = eye(m1_const + m1_row);

PI_sym = [
    eye(nu)
    (I-Delta*D)\Delta*C 
    ];

pLFR_sym = [A B]*PI_sym;

if ~isnumeric(pLFR_sym)
    PI_sym = simplify(PI_sym);
    pLFR_sym = simplify(pLFR_sym);
end

end

function [A,B,C,D,Delta] = data(matlfr)

[D,C,B,A,blk] = lfrdata(matlfr);

s = numel(blk.names);
c = cell(s, 1);
for k = 1:s
    c{k} = sym(blk.names{k}) * eye(blk.desc(1:2,k)');
end
Delta = blkdiag(c{:});

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function test1_lfr_with_constant_block
%%

lfrs a b c 


A = [
    1/a + b*c , 1/(b+c+5)
    12+a*c^2  , a*b*c
    ];

tic
lfr2sym(A)
toc

end

function test1_lfr_without_constant_block
%%

lfrs a b c 


A = [
    a*b + b*c , 1/(b+c+5)
    12+a*c^2  , a*b*c
    ];

tic
lfr2sym(A)
toc

end