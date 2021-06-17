function [ret] = pcz_persist_print(persist, printargs, varargin)
%% Script pcz_persist_print
%
%  file:   pcz_persist_print.m
%  author: Peter Polcz <ppolcz@gmail.com>
%
%  Created on 2017.08.25. Friday, 13:07:38
%
%%

[h,str,varargin] = pcz_persist_handle_or_name(varargin{:});

figname = persist.timefl('fig', str, varargin{:});

print(h,figname,printargs{:});

pcz_dispFunction('File saved to:\n%s.%s', figname, printargs{1}(3:end));

fnimage = pcz_resolvePath(figname);
fnscript = pcz_resolvePath(matlab.desktop.editor.getActive().Filename);

command = sprintf('krusader --left %s --right %s', fnscript.dir, fnimage.dir);
clipboard('copy', command);

end