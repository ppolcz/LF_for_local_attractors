function [C_num, T_num, C_den, T_den] = pcoeffs(A,vars)
%% 
%  
%  file:   pcoeffs.m
%  author: Polcz Péter <ppolcz@gmail.com> 
%  
%  Created on 2016.05.25. Wednesday, 11:12:53
%  Modified on 2019. October 31. (2019a)
%

if nargin < 2
    vars = symvar(A);
end

[num,den] = numden(A);

C_num = cell(size(A));
T_num = cell(size(A));

C_den = cell(size(A));
T_den = cell(size(A));

for i = 1:numel(A)
    [c,t] = coeffs(num(i),vars);
    C_num{i} = c;
    T_num{i} = t;
    
    [c,t] = coeffs(den(i),vars);
    C_den{i} = c;
    T_den{i} = t;
    
    C_num{i} = C_num{i} / C_den{i}(1);
    C_den{i} = C_den{i} / C_den{i}(1);
end


end