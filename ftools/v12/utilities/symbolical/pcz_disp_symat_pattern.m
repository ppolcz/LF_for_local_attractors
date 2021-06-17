function [ret] = pcz_disp_symat_pattern(A)
%%
%  File: pcz_disp_symat_pattern.m
%  Directory: 7_ftools/ftools/v12.01/utilities/symbolical
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2021. January 15. (2020b)
%


[n,m] = size(A);

pattern = repmat(' ',n,2*m+1);

for i = 1:n
    
    for j = 1:m
   
        if pcz_symzero1(A(i,j))
            pattern(i,2*j) = '0';
        else
            pattern(i,2*j) = '#';
        end
    end
    
end

disp(['Sparsity pattern for ' inputname(1) ':'])
disp ' '
for i = 1:n
    disp(['    [' pattern(i,:) ']'])
end
disp ' '

end