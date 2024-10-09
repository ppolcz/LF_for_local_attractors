function [ret] = pcz_num2str_v2_ginput_text(args)
arguments
    args.format = '%g';
    args.del1 = ',';
    args.del2 = '; ';
    args.del2end = '';
    args.pref = '';
    args.beginning = 'text(';
    args.ending = ',"");';
    args.label = '';
    args.round = 4;
    args.name = '[noinputname]';
    args.inputname = '{inputname}';    
end
%%
%
%  file:   pcz_num2str.m
%  author: Polcz Péter <ppolcz@gmail.com>
%
%  Created on 2017.02.12. Sunday, 20:49:45
%
% Examples:
% 
%  pcz_num2str(A,B,2*C,A+B,'format','%f', 'del1', ' , ', 'del2', ' ; ',...
%      'pref', '', 'beginning', ' [ ', 'ending', ' ] ', 'round', 5,...
%      'label', '{inputname} = \n', 'name', '2*C')
% 
%  pcz_num2str(A,b,'label','{inputname} = ', 'pref', '    ', 'del2','\n','beg', '[\n')
%
%  A = [
%      8 , 1 , 6
%      3 , 5 , 7
%      4 , 9 , 2 ]
%  b = [
%      0.8147 , 0.6324 , 0.9575 , 0.9572 , 0.4218
%      0.9058 , 0.0975 , 0.9649 , 0.4854 , 0.9157
%      0.127 , 0.2785 , 0.1576 , 0.8003 , 0.7922
%      0.9134 , 0.5469 , 0.9706 , 0.1419 , 0.9595 ]

qty = num2cell(pginput,2);

if ~iscell(qty)
    qty = {qty};
end

if nargout > 0
    ret = cell(numel(qty),1);
end

for i = 1:numel(qty)
    args.var = qty{i};
    if args.round > 0
        args.var = round(args.var,args.round);
    end

    if isempty(args.var)
        str = '[]';
    else
        str = pcz_num2str_fixed(args.var,args.format, args.del1, args.del2, args.pref, ...
            args.beginning, args.ending);
    end
        
    if nargout > 0
        ret{i} = str;
    else
        disp(str)
    end

end

if nargout > 0
    if numel(ret) == 1
        ret = ret{1};
    end
end

end