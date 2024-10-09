%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First add manually:                                                     %
% - MOSEK Solver (https://www.mosek.com/downloads/list/9/)                %
% - SMAC LFR Toolbox v2.1 (https://w3.onera.fr/smac/lfrt)                 %
% - SMAC GSS Library v2.1 (https://w3.onera.fr/smac/gss_download)         %
% to the Matlab path.                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ROOT = fileparts(mfilename('fullpath'));
cd(ROOT)

rmpath(genpath(ROOT))
