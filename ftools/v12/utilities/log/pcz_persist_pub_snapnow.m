function pcz_persist_pub_snapnow(persist, varargin)
%% Script pcz_persist_pub_snapnow
%
%  file:   pcz_persist_pub_snapnow.m
%  author: Peter Polcz <ppolcz@gmail.com>
%
%  Created on 2017.08.25. Friday, 12:25:39
%
%%

[h,str,varargin] = handle_or_name(varargin{:});

figname = [ persist.pub_absdir '/' str '.png' ];
html_figname = [ persist.pub_dirname '/' str '.png' ];
title = '';
if numel(varargin) > 0
    title = varargin{1};
end

drawnow
snapnow
pcz_savefig(h,figname,'fig',0,'export',1);
pause(0.1)

system(['mkdir -p ' persist.pub_absdir_thumb]);
thumb_cmd = [
    'mogrify -resize x170 -background white -gravity center ' ...
    '-extent x170 -format png -quality 75 -path ' persist.pub_absdir_thumb ' ' ...
    figname
    ];
system(thumb_cmd)

fieldname = sprintf('snapshot_html_%s', datestr(now,'HH_MM_SS_FFF'));
persist.(fieldname) = [
    '<?php thumb_gallery("' html_figname '", "' persist.file.fn ' snapshot", "' title '", "files/scripts"); ?>'
    ];
clipboard('copy', persist.(fieldname));

end