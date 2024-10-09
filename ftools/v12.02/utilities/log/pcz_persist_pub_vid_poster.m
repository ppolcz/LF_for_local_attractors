function pcz_persist_pub_vid_poster(persist, varargin)
%% Script pcz_persist_pub_vid_poster
%
%  file:   pcz_persist_pub_vid_poster.m
%  author: Peter Polcz <ppolcz@gmail.com>
%
%  Created on 2017.08.25. Friday, 12:24:36
%
%%

[h,str,~] = pcz_persist_handle_or_name(varargin{:});

persist.pub_vidname = str;
persist.runID

persist.pub_vid_poster_filename = sprintf('%s_poster.png', str);
figname = [ persist.pub_absdir '/' persist.pub_vid_poster_filename ];

if ~exist(persist.pub_absdir,'dir')
    mkdir(persist.pub_absdir)
end

drawnow
snapnow
pcz_savefig(h,figname,'fig',0,'export',1);
pause(0.1)

end