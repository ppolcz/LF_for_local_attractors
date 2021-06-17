function [ret] = sed_platex_from_latex
%% 
%  
%  file:   sed_platex_from_latex.m
%  author: Polcz Péter <ppolcz@gmail.com> 
%  
%  Created on 2016.02.20. Saturday, 01:49:07
%

ret = {
    % erase {ccc}
    %'\{c*\}' ''
    
    % replace \left(\begin{array} by \begin{bmatrix}
    %'\\\\left.\\\\begin\{array\}' '\\\\begin\{bmatrix\}'
    %'\\\\end\{array\}\\\\right.' '\\\\end\{bmatrix\}'
    
    % replace \left(\begin{array} by \begin{array}
    '\\\\left.\\\\begin\{array\}' '\\\\begin\{array\}'
    '\\\\end\{array\}\\\\right.' '\\\\end\{array\}'
    
    % replace \left \right with parathesis or brackets
    '\\\\(left)([^\\\\])' '\\\\qty\2'
    '\\\\(right)([^\\\\])' '\2'
    
    % erase the \mathrm{..} macros
    '\\\\mathrm\{([a-zA-Z0-9]*)\}' '\1'
    
    % erase the \, \! placeholder macros
    '\\\\,' ''
    '\\\\!' ''
    
    % in case of [a,b,c \\ c,d,e] --> .. a & b & c \\ c & d & e ..
    '([^\\\\]),' '\1 \& '
    
    % put before '\\' a space: ' \\'
    '(\\\\\\\\)' ' \1'

    % {x_2}^2 --> x_2^2
    '\{([^\{\}]*)\}(\^([a-zA-Z0-9]|\{[^\}]*\}))' '\1\2'
    
    % https://tio.run/#sed-gnu
    
    % {x_{2}}^2 --> x_{2}^2    
    '\{([a-zA-Z0-9]+)_\{([a-zA-Z0-9]+)\}\}\^([a-zA-Z0-9]+)' '\1_{\2}^\3'

    % \frac{\partial}{\partial x} f(x) --> \frac{\partial f(x)}{\partial x}
    '(\\\\frac\{\\\\partial)(\}\{\\\\partial [^\}]*\}) ([^ ]*)' '\1 \3\2'
    
    % remove extra space between sign operator and a variables
    '([\&\\\\]) *- +([a-zA-Z0-9])' '\1 -\2' % after & \
    '([\(]) *- +([a-zA-Z0-9])' '\1-\2' % after (

    };

end