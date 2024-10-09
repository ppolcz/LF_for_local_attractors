function [ret] = proot(scriptname, reldir, varargin)
%%
%
%  file:   proot.m
%  author: Polcz Péter <ppolcz@gmail.com>
%
%  Created on 2016.01.19. Tuesday, 22:14:05
%
%%

msgid = 'MATLAB:MKDIR:DirectoryExists';

ROOT = getenv('ROOT');
assert(~isempty(ROOT) && exist(ROOT,'dir'), ...
    'Root directory is not assigned well, first of all, run ./startup.m')

if nargin > 0 && ischar(scriptname)

    % This is a very ugly implementation
    fname = strrep(scriptname, [ ROOT '/' ],'');
    if ~isempty(fname) && strcmp(fname(end-1:end),'.m')
        fname = fname(1:end-2);
    end

    if nargin > 1 && iscell(reldir)
        midpath = reldir{1};

        if midpath(end) == '/'
            midpath = midpath(1:end-1);
        end

        tags = strjoin(varargin, '.');

        if isempty(tags)
            fname = [ midpath '/' fname '/' strrep(reldir{2},'/','.') ];
        else
            fname = [ midpath '/' fname '@' tags '/' strrep(reldir{2},'/','.') ];
        end

        % if isempty(tags)
        %     fname = [ midpath '/' strrep(fname,'/','.') '/' strrep(reldir{2},'/','.') ];
        % else
        %     fname = [ midpath '/' strrep(fname,'/','.') '@' tags '/' strrep(reldir{2},'/','.') ];
        % end

        ret = [ROOT '/' fname];

        warning('off', msgid);
        mkdir(ret)
        warning('on', msgid);
    end
else
    ret = ROOT;
end

end