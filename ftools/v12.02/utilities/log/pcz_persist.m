function [persist] = pcz_persist(fname, pcaller, varargin)
%% 
%  
%  file:   pcz_persist.m
%  author: Polcz Péter <ppolcz@gmail.com> 
%  
%  Created on 2016.05.16. Monday, 15:30:30
%  Reviewed on 2017. August 25. [decentralization, fonts, print functions]
%  
%  Examples:
%  
%   %fname: full path of the actual file
%   pcz_cmd_fname('fname');
%   persist = pcz_persist(fname);
%   persist.backup();
% 
%   saver = persist.savefig(figh, 'fig_nr%d', 4)
%   saver({'A', A}, {'B', B, 2});
%  
%   pcz_save(persist.simple('mat',format,args...), 'variable1', 'var2')
% 
% Publish video:

%% Initialization

if nargin == 0 || isempty(fname) || (pversion >= 2016 && contains(fname,'/tmp/Editor/'))
    try
        fname = matlab.desktop.editor.getActive().Filename;
    catch
        fname = '';
    end
end
persist = struct;

if nargin < 2
    pcaller = [];
end
    
assignin('base', 'fname', fname);

%% Fonts

% TeX Gyre Termes Math

persist.font.axis08 = {'FontSize',8,'FontName','TeX Gyre Schola Math'};
persist.font.axis10 = {'FontSize',10,'FontName','TeX Gyre Schola Math'};
persist.font.axis12 = {'FontSize',12,'FontName','TeX Gyre Schola Math'};
persist.font.axis14 = {'FontSize',14,'FontName','TeX Gyre Schola Math'};
persist.font.axis18 = {'FontSize',18,'FontName','TeX Gyre Schola Math'};
persist.font.latex18c = {'interpreter','latex','FontSize',18,'HorizontalAlignment','center'};
persist.font.latex18 = {'interpreter','latex','FontSize',18};
persist.font.latex = {'interpreter','latex'};

%%

% global variables
G = pglobals;
tags = varargin;

if nargin > 0

    % TIME-INDEPENDENT DIRECTORY DEFINITIONS
    
    % cache directory for the mat files
    persist.bin = proot(fname, G.RELPATH_MAT_FILES, tags{:});
    persist.mat = persist.bin;
    
    % cache directory for the figures
    persist.fig = proot(fname, G.RELPATH_FIGURES, tags{:});
    
    % cache directory for the figures
    persist.txt = proot(fname, G.RELPATH_TEXT_FILES, tags{:});

    % filenames
    persist.file = pcz_resolvePath(fname);
    
    % TIME-DEPENDENT DEFINITIONS

    % reset or inherint from caller
    if isempty(pcaller) || ~isfield(pcaller,'file') || ~strcmp(persist.file.path, pcaller.file.path)
        reset
        pcz_persist_display_status(get_persist,'initialized')
    else
        inherit(pcaller)
        pcz_persist_display_status(get_persist,'reused (inherited)')
    end
else
    reset
    pcz_persist_display_status(get_persist,'initialized (file independent)')
end

if nargin > 0
    
    % backup script filename
    persist.bfn = @backup_file_path;
    
    % backup command, which save the entire script file
    persist.backup = @backup;
    persist.reset = @reset;
    
    % filename generator function handles
    persist.simple = @simple;    
    persist.timefl = @timefl;    
    persist.savefig = @savefig;
    persist.saveall = @saveall;
    
    % Using Matlab's `print` function
    persist.png = @(varargin) pcz_persist_print(get_persist,{'-dpng'},varargin{:});
    persist.pdf = @(varargin) pcz_persist_print(get_persist,{'-dpdf'},varargin{:});
    persist.print = @(str,varargin) pcz_print(timefl('fig', str),varargin{:});

    persist.pub_vid_poster = @(varargin) pcz_persist_pub_vid_poster(get_persist, varargin{:});
    persist.pub_vid_write = @(varargin) pcz_persist_pub_vid_write(get_persist, varargin{:});
    persist.pub_snapnow = @(varargin) pcz_persist_pub_snapnow(get_persist, varargin{:});
end

% directory for temporary files
persist.tmp = [ proot G.TMP_PATH ];
if ~exist(persist.tmp,'dir')
    mkdir(persist.tmp);
end

% registry
persist.reg = [ proot G.REG_PATH ];
if ~exist(persist.reg,'dir')
    mkdir(persist.reg);
end

% registry
persist.session = [ proot G.SESSION_PATH ];
if ~exist(persist.session,'dir')
    mkdir(persist.session);
end

% possible registry names
names = cellfun(@(s) {s.name}, num2cell(dir(persist.reg)));
persist.reg_names = names(0 == cellfun(@isempty,regexp(names,'\.mat$')));

if isfield(persist,'file')
    persist.pub_html = '';
    persist.pub_dirname = sprintf('%s_%s', persist.file.bname, persist.stamp);
    persist.pub_reldir = sprintf('%s/%s', G.RELPATH_PUB,persist.pub_dirname);
    persist.pub_absdir = [ G.WWW '/' persist.pub_reldir ];
    persist.pub_absdir_thumb = [ persist.pub_absdir '/thumb' ];
    persist.pub_vidname = '';
    persist.pub_vid_poster_filename = '';
    persist.pub_vid_webm_filename = '';
    
    persist.pub_backup_script_path = ...
        sprintf('%s/%s', persist.pub_absdir, persist.file.fn);
    persist.pub_backup_script_relpath = ...
        sprintf('%s/%s', persist.pub_reldir, persist.file.fn);
end

%% Inner methods

    function reset
        % regenerate timestamp
        persist.fdate = pcz_fancyDate();
        persist.date = [ datestr(now, 'yy_mm_dd_') 'Time' datestr(now, 'HHMMSS')];
        
        if isfield(persist,'file')
            persist.runID = pcz_runID(persist.file.bname,'');
        else
            persist.runID = pcz_runID('','');
        end
        persist.stamp = [ persist.date '_runID' num2str(persist.runID) ];
    end

    function inherit(p)
        persist.fdate = p.fdate;
        persist.date = p.date;
        persist.runID = p.runID;
        persist.stamp = p.stamp;
    end

    function backup
        copyfile(persist.file.path, persist.bfn());
        pcz_dispFunction('Script `%s` backuped', persist.file.bname)
    end

    function fn = simple(type, str, varargin)
        fn = [ persist.(type) '/' sprintf(str, varargin{:}) ];
    end

    function fn = timefl(type, str, varargin)
        inputargs = [str varargin];
        fn = [ persist.(type) '/' persist.stamp '_' sprintf(inputargs{:})];
    end

    function saveall(str, varargin)
        filename = timefl('mat', str, varargin{:});
        evalin('caller', sprintf('save(''%s'')', filename));
        backup;
    end

    function [varargout] = savefig(varargin)
        
        [h,str,varargin] = pcz_persist_handle_or_name(varargin{:});
        
        figname = timefl('fig', str, varargin{:});
        fprintf('Figure saved to:\n%s\n', figname);
        pcz_savefig(h, figname);
        backup;
        
        if nargout > 0
            varargout{1} = @txtsaver;
        end
        
        function txtsaver(varargin)
            filename = [ figname '-data.txt' ];
            pcz_txtsaver(filename, varargin{:});
        end
    end

    function bfn = backup_file_path
        bfn = sprintf('%s/%s_%s%s', persist.txt, persist.file.bname, persist.stamp, persist.file.ext);
    end

    % FIGYELEM, HA A BASE WORKSPACE-BEN MEG NINCS PERSIST, EZ NEM A
    % TENYLEGES PERSIST-ET ADJA VISSZA
    function p = get_persist
        p = persist;
        % fprintf('pcz_persist.m::get_persist, runID = %d\n', p.runID);
    end

end
