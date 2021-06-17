function [ret] = pcz_persist_pub_vid_write(persist, frames, caption)
%% Script pcz_persist_pub_vid_write
%
%  file:   pcz_persist_pub_vid_write.m
%  author: Peter Polcz <ppolcz@gmail.com>
%
%  Created on 2017.08.25. Friday, 12:23:28
%
%%

if nargin <= 2
    caption = '';
end

if isempty(persist.pub_vidname)
    persist.pub_vidname = sprintf('vid%s', pcz_runID(caption,'pub_vid_write'));
end

persist.pub_vid_webm_filename = [ persist.pub_vidname '.webm' ];

avipath = [ persist.pub_absdir '/' persist.pub_vidname '.avi'];
webmpath = [ persist.pub_absdir '/' persist.pub_vid_webm_filename ];

if ~exist(persist.pub_absdir,'dir')
    mkdir(persist.pub_absdir)
end

v = VideoWriter(avipath);
open(v);
writeVideo(v,frames);
close(v);

cmd = sprintf('avconv -i %s -vf scale=1298:-1 -acodec copy -threads 4 %s',...
    avipath, webmpath);
try
    [status,results] = system(cmd);
catch ex
    error(getReport(ex))
end

nl = newline;

html = [ ...
    '<video width="620" style="max-width:100%" poster="' persist.pub_vid_poster_filename '" controls>' nl ...
    '    <source src="' persist.pub_vid_webm_filename '" type=''video/webm; codecs="vp8.0, vorbis"''>' nl ...
    '    Your browser does not support the video tag.' nl ...
    '</video>' ];
clipboard('copy', html);

if nargout > 0
    ret = html;
end

fprintf('|%s|\n',persist.pub_vid_webm_filename)

end