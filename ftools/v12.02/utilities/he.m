function [He,O,I] = he
%%
%  File: he.m
%  Directory: 7_ftools/ftools/v12/utilities
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2020. June 12. (2019b)
%

He = He_Class;
O = @(varargin) zeros(varargin{:});
I = @(varargin) eye(varargin{:});
end
