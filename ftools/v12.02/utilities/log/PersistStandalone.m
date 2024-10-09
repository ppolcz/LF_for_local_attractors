classdef PersistStandalone
    %%
    %
    %  file:   PersistStandalone.m
    %  author: Polcz Péter <ppolcz@gmail.com>
    %
    %  Created on 2018.02.11. Sunday, 11:17:40
    %

    %%
    properties (Constant = true)
        font_axis08 = {'FontSize',8,'FontName','TeX Gyre Schola Math'};
        font_axis10 = {'FontSize',10,'FontName','TeX Gyre Schola Math'};
        font_axis12 = {'FontSize',12,'FontName','TeX Gyre Schola Math'};
        font_axis14 = {'FontSize',14,'FontName','TeX Gyre Schola Math'};
        font_axis18 = {'FontSize',18,'FontName','TeX Gyre Schola Math'};
        font_axis22 = {'FontSize',22,'FontName','TeX Gyre Schola Math'};
        font_axis26 = {'FontSize',26,'FontName','TeX Gyre Schola Math'};
        font_latex18c = {'interpreter','latex','FontSize',18,'HorizontalAlignment','center','FontName','TeX Gyre Schola Math','Color',[0 0 0]};
        font_latex8  = {'interpreter','latex','FontSize', 8,'FontName','TeX Gyre Schola Math','Color',[0 0 0]};
        font_latex10 = {'interpreter','latex','FontSize',10,'FontName','TeX Gyre Schola Math','Color',[0 0 0]};
        font_latex12 = {'interpreter','latex','FontSize',12,'FontName','TeX Gyre Schola Math','Color',[0 0 0]};
        font_latex14 = {'interpreter','latex','FontSize',14,'FontName','TeX Gyre Schola Math','Color',[0 0 0]};
        font_latex16 = {'interpreter','latex','FontSize',16,'FontName','TeX Gyre Schola Math','Color',[0 0 0]};
        font_latex18 = {'interpreter','latex','FontSize',18,'FontName','TeX Gyre Schola Math','Color',[0 0 0]};
        font_latex22 = {'interpreter','latex','FontSize',22,'FontName','TeX Gyre Schola Math','Color',[0 0 0]};
        font_latex26 = {'interpreter','latex','FontSize',26,'FontName','TeX Gyre Schola Math','Color',[0 0 0]};
        font_latex = {'interpreter','latex'};
    end

    methods(Static)
        function ret = latexify_axis(ax, ax_fontsize)
            if isnumeric(ax)
                ax_fontsize = ax;
                ax = gca;
            end
            
            set(ax.Parent, 'Color', [1 1 1])
            
            set(ax, 'FontName','TeX Gyre Schola Math',...
                'GridColor', [0.1 0.1 0.1], 'MinorGridColor', [0.1 0.1 0.1],...
                'XColor', [0 0 0], 'YColor', [0 0 0], 'ZColor', [0 0 0]);
                        
            set(ax,'TickLabelInterpreter', 'latex');
            
            xl = get(ax,'XLabel');
            xlFontSize = get(xl,'FontSize');
            xAX = get(gca,'XAxis');
            set(xAX,'FontSize', ax_fontsize)
            set(xl, 'FontSize', xlFontSize);

            yl = get(ax,'YLabel');
            ylFontSize = get(yl,'FontSize');
            yAX = get(gca,'YAxis');
            set(yAX,'FontSize', ax_fontsize)
            set(yl, 'FontSize', ylFontSize);

            zl = get(ax,'ZLabel');
            zlFontSize = get(zl,'FontSize');
            zAX = get(gca,'ZAxis');
            set(zAX,'FontSize', ax_fontsize)
            set(zl, 'FontSize', zlFontSize);
            
            if nargout > 0
                ret = ax;
            end
        end
        
        function ret = latexify_colorbar(cbar, fontsize)
            set(cbar, 'FontSize', fontsize, 'FontName','TeX Gyre Schola Math',...
                'TickLabelInterpreter', 'latex');
            cbar.Label.Color = [0 0 0];
            cbar.Label.Interpreter = 'latex';

            if nargout > 0
                ret = cbar;
            end
        end
        
        function ret_ = latexified_labels(ax, fontsize, xl, yl, zl)
            
            ret = cell(1,nargin - 2);
            if nargin > 2
                ret{1} = xlabel(ax,xl,'interpreter','latex','FontSize', fontsize,'FontName','TeX Gyre Schola Math','Color',[0 0 0]);
                
            end
                
            if nargin > 3
                ret{2} = ylabel(ax,yl,'interpreter','latex','FontSize', fontsize,'FontName','TeX Gyre Schola Math','Color',[0 0 0]);
            end
            
            if nargin > 4
                ret{3} = zlabel(ax,zl,'interpreter','latex','FontSize', fontsize,'FontName','TeX Gyre Schola Math','Color',[0 0 0]);
            end
            
            % ret{1}.VerticalAlignment = 'middle';
            % ret{2}.VerticalAlignment = 'middle';

            if nargout > 0
                ret_ = ret;
            end
        end
        
        function path = generate_path(fname,dir,varargin)
            
            abspath = dir{1};

            tags_period = strjoin(varargin, '.');

            if isempty(tags_period)
                path = [ abspath '/' fname '/' strrep(dir{2},'/','.') ];
            else
                path = [ abspath '/' fname '@' tags_period '/' strrep(dir{2},'/','.') ];
            end

            path = cleanpath(path);
            
            msgid = 'MATLAB:MKDIR:DirectoryExists';
            warning('off', msgid);
            mkdir(path)
            warning('on', msgid);
    
        end

    end
    
    properties (GetAccess = public, SetAccess = protected)
        G;

        fname, fname_active;

        fdate, date, runID, runID_dates, stamp;
    end

    properties (GetAccess = public, SetAccess = public)
    end

    methods (Access = public)

        %%%
        % Constructor
        %
        function p = PersistStandalone()
            try
                p.fname_active = matlab.desktop.editor.getActive().Filename;
            catch
                p.fname_active = '';
            end

            p.fname = cleanpath(p.fname_active);

            if endsWith(p.fname,'.m')
                p.fname = p.fname(1:end-2);
            end
            
            assignin('base', 'fname', p.fname);
            assignin('base', 'fname_active', p.fname_active);

            p.G = pglobals;
            
            % p.fname is relative to the root path
            p.fname = strrep(p.fname,cleanpath([ p.G.ROOT '/' ]),'');
            
        end


        %%%
        % Regenerate timestamp
        %
        function p = reset(p)
            p.fdate = pcz_fancyDate();
            p.date = [ datestr(now, 'yyyy-mm-dd_') datestr(now, 'HH-MM-SS')];

            p = p.generate_runID;

            p.stamp = [ p.date '_i' num2str(p.runID) ];
        end


        %%%
        % Generate runID (overrided by Persist)
        %
        function p = generate_runID(p)
            [p.runID,p.runID_dates] = pcz_runID('','');
            p.display_status;
        end


        %%%
        % Display status (overrided by Persist)
        %
        function display_status(p)
           pcz_dispFunction('Persistence initialized (file independent) [run ID: %d, %s]', ...
                p.runID, p.fdate);
        end
    end

    methods (Access = private)

    end



end
