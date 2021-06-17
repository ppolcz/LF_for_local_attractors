function [cmdout,status] = sed_apply(input, rules)
%% 
%  
%  file:   sed_apply.m
%  author: Polcz Péter <ppolcz@gmail.com> 
%  
%  Created on 2016.02.20. Saturday, 01:56:46
%

warning off MATLAB:lang:cannotClearExecutingFunction

% eg. sed_cmd = s/replace_this/to_this/g;..;s/kutyagumi/valami/g
rules_str = strjoin(cellfun(@(a,b) {['s/' a '/' b '/g']}, rules(:,1), rules(:,2))', ';');
[status,cmdout] = system(['printf ''%s'' ''' input ''' | sed -r "' rules_str '"']);

warning on MATLAB:lang:cannotClearExecutingFunction


end