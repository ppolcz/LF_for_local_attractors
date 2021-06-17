function [t,p,dp] = pcz_pdp_lim_smooth_v2(T, dt, p_lim, dp_lim, varargin)
%% pcz_pdp_lim_smooth_v2
%  
%  File: pcz_pdp_lim_smooth_v2.m
%  Directory: 2_demonstrations/lib/matlab
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. November 11. (2019a)
%

%%

% Nr. of parameters
np = size(p_lim,1);

pl = p_lim * [-1;1];
dpl = min(abs(dp_lim),[],2);

T_ = T * dpl ./ pl;

t = ( 0:dt:T )';

N = numel(t);

p = zeros(N,np);
dp = zeros(N,np);

for i = 1:np
    [ti,pi,dpi] = pcz_pdp_lim_smooth_normalized(T_(i),dt,varargin{:});
    
    p(:,i) = interp1(ti*T/T_(i),pi,t) * pl(i) + p_lim(i,1);
    dp(:,i) = interp1(ti*T/T_(i),dpi,t) * T_(i)/T * pl(i);
    
end

end

function test1
%%

p_lim = [
    -1 1
    0 1
    0 1
    ];

dp_lim = [
    -2 2
    -1 1
    -100 100
    ];

np = numel(p_lim)/2;
T = 10;
dt = 0.01;

[t,p,dp] = pcz_pdp_lim_smooth_v2(T,dt,p_lim,dp_lim,'eltoldp',false,'plot',false);

pl(t,p,dp,p_lim,dp_lim,np)

end

function pl(t,p,dp,p_lim,dp_lim,np)

    figure(1)
    for i = 1:np
        ax(1) = subplot(2,np,i); plot(t,p(:,i)), grid on
        title(sprintf('p%d',i))
        xlim(t([1,end]))
        ylim(p_lim(i,:))

        grid on
        ax(2) = subplot(2,np,np+i); plot(t,dp(:,i)), grid on
        title(sprintf('dp%d',i)) 
        xlim(t([1,end]))
        ylim(dp_lim(i,:))

        linkaxes(ax,'x')
    end

end