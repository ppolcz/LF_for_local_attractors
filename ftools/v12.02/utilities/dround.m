function [ret] = dround(x, prec)
%% 
%  
%  file:   pround.m
%  author: Polcz Péter <ppolcz@gmail.com> 
%  
%  Created on 2016.04.19. Tuesday, 14:43:23
%

% axiliary function handle for numeric simplification
mag = ceil(log10(abs(x)));         % 10^(mag-1) < x <= 10^mag
mul = 10.^(prec-mag);              % multiplier in the round
ret = round(x.*mul) ./ mul;        % round to 

end