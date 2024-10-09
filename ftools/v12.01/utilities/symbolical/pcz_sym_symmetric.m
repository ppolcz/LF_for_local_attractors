function [P,p,p_idx] = pcz_sym_symmetric(name,n,~)
%% pcz_sym_symmetric
%  
%  File: pcz_sym_symmetric.m
%  Directory: 2_demonstrations/lib/matlab
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. October 23.
%

%%

nr = n*(n+1)/2;

p = pcz_sym_indexed(name,nr).';

U = triu(ones(n));

p_idx = find(U == 1);

P = sym(U);

P(p_idx) = p;

P = P + triu(P,1).';

end

function self_check
%%
% TODO: erre nem mukodik

P = pcz_sym_symmetric('p',4)

end
