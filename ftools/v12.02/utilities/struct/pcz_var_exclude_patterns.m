function [exclude] = pcz_var_exclude_patterns(varargin)
%% 
%  
%  file:   pcz_var_exclude_patterns.m
%  author: Polcz Péter <ppolcz@gmail.com> 
%  
%  Created on 2017.02.14. Tuesday, 23:46:40
%

exclude = {
    '^TMP_*'
    '^props$'
    '^opts$'
    '^persist$'
    '^fname$'
    '.*_old_\d{2}_\d{2}_\d{2}'
    };

%     '^$'

end