function [G,globscode_] = pglobals
%% 
%  
%  file:   pcz_globals.m
%  author: Polcz Péter <ppolcz@gmail.com> 
%  
%  Created on 2016.01.19. Tuesday, 13:55:28
%
%% ROOT

ROOT = getenv('ROOT');
assert(~isempty(ROOT) && exist(ROOT,'dir'), ...
    'Root directory is not assigned well, first of all, run ./startup.m')

G = struct;
G.ROOT = ROOT;
G.WWW = getenv('WWW');
G.WWW_FILES = getenv('WWW_FILES');
G.WWW_HTMLS = getenv('WWW_HTMLS');
G.EXTERNAL_DIR = getenv('EXTERNAL_DIR');
G.EXTERNAL_LOG_DIR = getenv('EXTERNAL_LOG_DIR');
G.PYLIB = getenv('PYLIB');

G.STYLESHEET = [G.ROOT '/publish.xsl' ];

assert(~isempty(G.ROOT) && exist(G.ROOT,'dir'), sprintf('%s `%s` %s',...
    'Environmental variable',...
    'ROOT',...
    'should exists and it should point to a valid path!'))

%% 

% External loggind directories
G.RELPATH_MAT_FILES = { G.EXTERNAL_LOG_DIR 'binary' };
G.RELPATH_FIGURES = { G.EXTERNAL_LOG_DIR 'media' };
G.RELPATH_TEXT_FILES = { G.EXTERNAL_LOG_DIR 'text' };

% External registries
G.REG_PATH = getenv('EXTERNAL_LOG_DIR');
G.TMP_PATH = getenv('EXTERNAL_TMP_DIR');
G.SESSION_PATH = getenv('EXTERNAL_SESSION_DIR');

%% The A registry

G.REG_A = getenv('A_MAT_PATH');
G.SAVE_A = @(var) save_variable(G.REG_A,inputname(1),var);

%% publish

G.RELPATH_PUB = G.WWW_FILES;
G.VIEW_SCRIPTS = G.WWW_HTMLS;
G.PUB = G.WWW_FILES;

G.SANITIZER = getenv('MATLAB_HTML_SANITIZER');

G.PUB_SUBL = strcmp('true',getenv('PUB_SUBL'));
G.PUB_WEBVIEW = strcmp('true',getenv('PUB_WEBVIEW'));
G.PUB_EDIT = strcmp('true',getenv('PUB_EDIT'));
 
%%

globscode = sprintf('G_reset\n');

if nargout > 1
    globscode_ = globscode;
end

if nargout == 0
    clipboard('copy',globscode);
end

end

function save_variable(mat,varname,var)
    s.(varname) = var;
    save(mat,'-struct','s');
end
