function [c_ret,t_ret] = pcz_coeffs(expression, variables)
%% 
%  
%  file:   pcz_coeffs.m
%  author: Polcz Péter <ppolcz@gmail.com> 
%  
%  Created on 2016.03.06. Sunday, 22:06:05
%

if nargin == 2

    warning 'This branch of the function is a little bit buggie!'
    % Try these:
    % syms a b c d real
    % A = [a^2+b]
    % 1) v = [a b c d a^2]' , 2) v = [a b c d]' , 3) v = [a b c d]
    % pcz_coeffs(A, v)
    
    n = numel(variables);

    [q,m] = size(expression);
    t = 1+n;

    expr_coeffs = sym(zeros(q,m*t));
    expr_0 = subs(expression,variables,variables*0);

    % [1] generate the coefficients matrix
    E = [1 ; zeros(n-1,1) ];
    for i = 1:n
        expr_coeffs(:,(1:t:m*t)+i) = subs(expression, variables, E) - expr_0;
        E = circshift(E,1);
    end

    expr_coeffs(:,1:t:m*t) = expr_0;

    if nargout > 0
        c_ret = expr_coeffs;
    end
        
elseif nargin == 1

    if ~isscalar(expression)
        error 'Input argument `expression` should be a scalar variable'
    end
    
    if isempty(symvar(expression))
        try
            tmp = double(expression); %#ok
        catch
            assert(false, 'serious problem detected')
        end
        [c,t] = deal(expression, sym(1));
    else
        
        [c,t] = coeffs(expression);
        
    end
    
    if nargout > 0
        c_ret = c;
    end
    
    if nargout > 1
        t_ret = t;
    end
end

end