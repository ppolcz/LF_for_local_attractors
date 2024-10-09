classdef Persist < PersistStandalone
    %%
    %
    %  file:   Persist.m
    %  author: Polcz Péter <ppolcz@gmail.com>
    %
    %  Created on 2018.02.11. Sunday, 00:03:09
    %

    %%
    properties (Constant = true)
    end

    properties (GetAccess = private, SetAccess = private)
    end

    properties (GetAccess = public, SetAccess = private)
        tags;
        pcaller = [];

        bin, mat, fig, txt, log, file;
        
        % Backup file's name
        bfn;
        
        % Diary file's name
        dfn;
        
        startTime;

        b_inherit;

        logfile = '';
        
        pub_html;
        pub_dirname;
        pub_reldir;
        pub_absdir;
        pub_absdir_thumb;
        pub_vidname;
        pub_vid_poster_filename;
        pub_vid_webm_filename;
        pub_backup_script_path;
        pub_backup_script_relpath
    end

    properties (GetAccess = public, SetAccess = public)
        simname = 'results';        
    end

    methods (Access = public)

        %%%
        % Constructor
        %
        function p = Persist(fname, pcaller, varargin)
            p@PersistStandalone()

            p.b_inherit = nargin >= 2;

            if ~contains(fname,'/tmp/Editor/')
                p.fname = strrep(cleanpath(fname),cleanpath([ p.G.ROOT '/' ]),'');
            end
                
            p.tags = varargin;
            assignin('base', 'fname', p.fname);

            if p.b_inherit
                p.pcaller = pcaller;
            end

            p = p.timeindep_dir_definitions;
            p = p.reset_or_inherint_from_caller;
            p = p.timedep_path_definitions;
            p = p.publish_path_definitions;
            
            % Start diary (output logging)
            diary(p.dfn);
            pcz_info('Output logging (with `diary''): %s', p.dfn);
            
            p.startTime = pcz_dispFunctionName;
        end


        %%%
        % Generate script and tag dependent folder names. Date values are
        % not used here.
        %
        function p = timeindep_dir_definitions(p)
            % cache directory for the mat files
            p.bin = Persist.generate_path(p.fname, p.G.RELPATH_MAT_FILES, p.tags{:});
            p.mat = p.bin;

            % cache directory for the figures
            p.fig = Persist.generate_path(p.fname, p.G.RELPATH_FIGURES, p.tags{:});

            % cache directory for the figures
            p.txt = Persist.generate_path(p.fname, p.G.RELPATH_TEXT_FILES, p.tags{:});

            % cache directory for the figures
            p.log = Persist.generate_path(p.fname, p.G.RELPATH_TEXT_FILES, p.tags{:});

            % filenames
            p.file = pcz_resolvePath(p.fname);
        end

        
        function p = timedep_path_definitions(p)
            p.bfn = sprintf('%s/%s_%s%s', p.txt, p.file.bname, p.stamp, p.file.ext);
            p.dfn = sprintf('%s/%s_%s%s', p.txt, p.file.bname, p.stamp, '_output.txt');
        end


        function p = publish_path_definitions(p)
            p.pub_html = '';
            p.pub_dirname = sprintf('%s_%s', p.file.bname, p.stamp);
            p.pub_reldir = sprintf('%s/%s', p.G.RELPATH_PUB,p.pub_dirname);
            p.pub_absdir = [ p.G.WWW '/' p.pub_reldir ];
            p.pub_absdir_thumb = [ p.pub_absdir '/thumb' ];
            p.pub_vidname = '';
            p.pub_vid_poster_filename = '';
            p.pub_vid_webm_filename = '';

            p.pub_backup_script_path = ...
                sprintf('%s/%s', p.pub_absdir, p.file.fn);
            p.pub_backup_script_relpath = ...
                sprintf('%s/%s', p.pub_reldir, p.file.fn);
        end

        %%%
        %
        %
        function p = reset_or_inherint_from_caller(p)
            % reset or inherint from caller
            if isempty(p.pcaller) || ~strcmp(p.file.path, p.pcaller.file.path)
                p = p.reset;
                p.display_status('initialized');
            else
                p = p.inherit;
                p.display_status('reused (inherited)');
            end
        end

        function p = inherit(p)
            p.fdate = p.pcaller.fdate;
            p.date = p.pcaller.date;
            p.runID = p.pcaller.runID;
            p.stamp = p.pcaller.stamp;
        end

        function p = generate_runID(p)
            [p.runID,p.runID_dates] = pcz_runID(p.file.bname,'');
        end

        function backup(p)
            if ~exist(p.file.path, 'file'), return, end
            copyfile(p.file.path, p.bfn);
            pcz_info('Script `%s` backed up.', p.file.bname)
        end

        function fn = file_simple(p, type, str, varargin)
            fn = [ p.(type) '/' sprintf(str, varargin{:}) ];
        end

        function fn = file_timefl(p, type, str, varargin)
            fn = [ p.(type) '/' p.stamp '_' sprintf(str, varargin{:})];
        end


        %%%
        % Save all workspace variables
        %
        function saveall(p, str, varargin)
            filename = p.file_timefl('mat', str, varargin{:});
            evalin('caller', sprintf('save(''%s'')', filename));
            p.backup;
        end


        %%%
        % Save certain variables.
        % Created on 2018.04.04. (április  4, szerda), 17:27
        %
        function save(p, filename, varargin)
            [dname,bname,~] = fileparts(filename);

            if ~isempty(dname)
                warning('relative path should be empty! Using only the basename!')
            end

            filename = p.file_timefl('mat', [ bname '.mat' ]);

            [dname,~,~] = fileparts(filename);
            
            if ~exist(dname, 'dir')
                system(sprintf('mkdir -p "%s"', dname))
                
                pcz_info('Directory %s not exists, created!', dname)
            end
                
            argnames = cell(1,nargin-1);
            argnames{1} = filename;
            
            for i = 3:nargin
                argnames{i-1} = inputname(i);
            end
            
            savedvars = join(argnames(2:end), ', ');
            
            argnames = cellfun(@(n) {['''' n '''']}, argnames);
            argnames = join(argnames,',');
            
            if exist(filename,'file')
                pcz_info('File already exists, append variables.')
                cmd = sprintf('save(%s,''-append'')', argnames{1});
            else
                cmd = sprintf('save(%s)', argnames{1});
            end
                
            evalin('caller', cmd);
            p.backup;
            
            pcz_info('Saved variables: %s.', savedvars{1});
            
            clipboard('copy',filename);
            
            pcz_info('Mat filed copied to clipboard.')
        end



        function [varargout] = savefig(p,varargin)

            [h,str,varargin] = pcz_persist_handle_or_name(varargin{:});

            figname = p.file_timefl('fig', str, varargin{:});
            fprintf('Figure saved to:\n%s\n', figname);
            pcz_savefig(h, figname);
            p.backup;

            if nargout > 0
                varargout{1} = @txtsaver;
            end

            function txtsaver(varargin)
                filename = [ figname '-data.txt' ];
                pcz_txtsaver(filename, varargin{:});
            end
        end


        %%%
        % Display status (PersistStandalon is overrided)
        %
        function display_status(p, msg)
            if p.b_inherit
                pcz_info('Persistence for `%s` %s [run ID: %s, %s]', ...
                    p.file.bname, msg, p.runID, p.fdate);
            else
                pcz_info('Persistence - empty %s [run ID: %s, %s]', ...
                    msg, p.runID, p.fdate);
            end
        end


        %%%
        % Examples: (you can try also `pcz_print`)
        %
        %   persist.print('-dpng', 'main_figure.png')
        %   persist.print({'-dpng' '-r100'}, 'main_figure.png')
        %
        function print(p, printargs, varargin)
            [h,str,varargin] = Persist.handle_or_name(varargin{:});

            figname = p.file_timefl('fig', str, varargin{:});
            
            [dir,bname,ext] = fileparts(figname);
            
            figname_trimmed = [dir '/' bname '-trimmed' ext];

            if ~iscell(printargs)
                printargs = {printargs};
            end
            
            figname, printargs{:}
            
            print(h,figname,printargs{:})

            pcz_dispFunction('File saved to:\n%s', figname);

            fnimage = pcz_resolvePath(figname);
            fnscript = pcz_resolvePath(matlab.desktop.editor.getActive().Filename);

            command = sprintf('krusader --left %s --right %s', fnscript.dir, fnimage.dir);
            clipboard('copy', command);
            
            if strcmp(ext,'.png') || strcmp(ext,'.jpg')
                pcz_convert_trim(figname,figname_trimmed,...
                    'Border',10,'BorderColor','White');
            end
        end

        function stoplog(p)
            TMP = p.startTime;
            pcz_dispFunctionEnd(TMP);
            
            diary off
            if exist(p.dfn,'file')
                pcz_output2log(p.dfn);
                pcz_dispFunction2(' ')
                pcz_info('Logfile formatted!');
            end
        end
        
    end
    
    methods (Static)
        
        function [h,str,varargin] = handle_or_name(h, str, varargin)
            if nargin < 2 || isempty(str)
                str = '';
            end

            if ischar(h) || iscell(h)
                varargin = [str varargin];
                str = h;
                h = gcf;
            end
        end
        
    end

    methods (Access = private)

        function obj = init_fonts(obj)
        end

    end



end
