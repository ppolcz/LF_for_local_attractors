function [ret] = pcz_latex_v7(A, varargin)
%% pcz_dispFunctionLatex
%  
%  File: pcz_dispFunctionLatex.m
%  Directory: 2_demonstrations/lib/matlab
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. May 30.
%

%%

s.label = inputname(1);
s = parsepropval(s,varargin{:});

str = '';

if ~isnumeric(A)
    str = pcz_latex_v6(A, 'disp_mode', 2, 'label', [ s.label ' = ' ], varargin{:});

    % 2021.03.08. (március  8, hétfő), 12:26
    % if size(A,2) > 10
    cccccc = repmat('c',[1 size(A,2)]);        
    str = strrep(str, ['\left(\begin{array}{' cccccc '}'], '\pmqty{');
    str = strrep(str, '\end{array}\right)', '}');
    % end
    
elseif isnumeric(A)
    tmp.A = A;
    str = pcz_num2str_latex(tmp.A,'beg',[s.label ' = \\spmqty{ '],'end',' }');    
end

disp(str)

end
