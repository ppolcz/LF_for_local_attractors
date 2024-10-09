function [ret] = pcz_cat_matrix(F)
%% pcz_cat_matrix
%  
%  File: pcz_cat_matrix.m
%  Directory: utilities/matrix_operations
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2020. May 14. (2019b)
%

%%

ret = cellfun(@(col) {vertcat(col{:})}, num2cell(F,1));
ret = [ ret{:} ];

end