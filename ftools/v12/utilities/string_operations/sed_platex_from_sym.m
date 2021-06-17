function [ret] = sed_platex_from_sym
%% 
%  
%  file:   set_platex_from_sym.m
%  author: Polcz Péter <ppolcz@gmail.com> 
%  
%  Created on 2016.02.20. Saturday, 01:49:35
%

ret = {
    % replace a1_{2} by a_{12}
    '([a-zA-Z])([0-9])_\{([0-9])\}' '\1_{\2\3}'
    
    % replace a1 by a_1 or s9 by s_9
    '([a-zA-Z])([0-9])' '\1_\2'
    
    % replace a12 by a_{12}
    '([a-zA-Z])([0-9][0-9]+)' '\1_{\2}'
    
    };

end