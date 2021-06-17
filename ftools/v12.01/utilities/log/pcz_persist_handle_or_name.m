function [h,str,varargin] = pcz_persist_handle_or_name(h, str, varargin)
%% Script pcz_persist_handle_or_name
%
%  file:   pcz_persist_handle_or_name.m
%  author: Peter Polcz <ppolcz@gmail.com>
%
%  Created on 2017.08.25. Friday, 12:55:55
%
%%

if nargin < 2 || isempty(str)
    str = '';
end

if ischar(h) || iscell(h)
    varargin = [str varargin];
    str = h;
    h = gcf;
end

end