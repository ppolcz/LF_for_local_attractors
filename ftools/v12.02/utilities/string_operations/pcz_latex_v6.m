function [ret] = pcz_latex_v6(A, varargin)
%% Script pcz_latex_v6
%  
%  File: pcz_latex_v6.m
%  Directory: projects/3_outsel/2018_01_10_LPV_inversion
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. January 13.
%
%%

o.label = '';
o.sed_user = {};
o.disp_mode = 2;
o.prec = 0;
o.vline = [];
o.hline = [];
o = parsepropval(o,varargin{:});

if ~isempty(o.hline), o.disp_mode = 2; end

label = sprintf(o.label);

disp_mode = o.disp_mode;
sed_user = o.sed_user;
prec = o.prec;

indent = '    ';

if iscell(disp_mode)
    sed_user = disp_mode;
    disp_mode = 0;
end

assert(isscalar(disp_mode), 'Arguments are wrong: display mode should be scalar')

sed = [
    sed_platex_from_latex
    sed_platex_from_sym
    ];

if disp_mode >= 1
    sed = [ sed ; {
        
        % newline+indent after \begin{..} or \begin{..}{..}, 
        % after \\, before \end{..}
        '(\\\\begin\{[^\}]*\}(\{[^\}]*\})*)\\s*' ['\1\n' indent]
        '(\\\\\\\\)' ['\1\n' indent]
        '(\\\\end\{)' '\\\\\\\\\n\1'
        
        % handle/erase extra space at the beggining of the lines
        [indent '*- '] [indent '-']
        [indent '*-'] [indent '-']
        [indent ' *'] indent
        } ];
    if disp_mode >= 3
        sed = [ sed ; {
            '(\\\\\\\\)' ['\n' indent '\1']
            '(&)' ['\n' indent '\1\n' indent]
            
            % remove extra newline before \end{..}
            [indent '\\\\\\\\\n(\\\\end)'] '\1'
            } ];
    end
end

sed = [ sed ; sed_user ];

try
    A_latex = latex(A);
catch
    try
        A_latex = pcz_matrix2latex_3(A, prec);
    catch ex
        getReport(ex)
    end
end

cccccc = repmat('c',[1 size(A,2)]);
coldef = cccccc;
if ~isempty(o.vline)
    bar = repmat(' ',[1 size(A,2)+1]);
    bar(o.vline+1) = '|';
    coldef = [ bar(1) reshape([coldef ; bar(2:end)], [1 size(A,2)*2]) ];
    coldef = strrep(coldef, ' ', '');
end 
A_latex = strrep(A_latex, cccccc, coldef);
   
% escaped = strrep(A_latex,'\', '\\');
% [cmdout,status] = sed_apply(escaped, sed);

[cmdout,status] = sed_apply(A_latex, sed);

nl_index = strfind(cmdout, newline);
for linenr = sort(o.hline,'desc')
    i = nl_index(linenr+1);
    cmdout = [cmdout(1:i-1) ' \hline ' cmdout(i:end) ];
end

nl_index = strfind(cmdout, newline);
if (disp_mode == 2) && (numel(nl_index)-1 == size(A,1))
    lines = cell(1,size(A,1));
    places = zeros(size(A,1),size(A,2));
    
    for i = 1:size(A,1)
        lines{i} = cmdout(nl_index(i)+1:nl_index(i+1)-1);
        places(i,1:end-1) = strfind(lines{i},'&');
        places(i,end) = strfind(lines{i}, '\\');
    end
    
    for j = 1:size(A,2)
        max_pos = max(places(:,j));
        for ii = 1:size(A,1)
            pos = places(ii,j);
            xspace = repmat(' ', [1 max_pos-pos]);
            lines{ii} = [lines{ii}(1:pos-1) xspace lines{ii}(pos:end)];
            places(ii, j+1:end) = places(ii, j+1:end) + max_pos-pos;
        end
    end
    
    cmdout = [ cmdout(1:nl_index(1)) strjoin(lines,newline) cmdout(nl_index(end):end) ];
end

if ~isscalar(A)
    cmdout = [ '\left(' cmdout '\right)' ];
end

if ~isempty(inputname(1)) && isempty(o.label)
    name = [inputname(1) ' = '];
else
    name = '';
end

cmdout = [label name cmdout];

if status
    fprintf(['\nreturn_value = ' num2str(status) '\n\n'])
end

if nargout > 0
    ret = cmdout;
else
    if disp_mode >= 1
        disp ' '
    end
    disp(cmdout)
end


end