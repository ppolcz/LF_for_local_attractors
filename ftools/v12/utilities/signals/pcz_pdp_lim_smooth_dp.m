function [t,p,dp] = pcz_pdp_lim_smooth_dp(T, dt, p_lim, dp_lim, varargin)
%% pcz_pdp_lim_smooth_v2
%  
%  File: pcz_pdp_lim_smooth_v2.m
%  Directory: 2_demonstrations/lib/matlab
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. November 09. (2019a)
%

TMP_bzqpLWoBcwMvZxPDPVtb = pcz_dispFunctionName;

%%

opts.nrint = 4;
opts.maxnrit = 10;
opts.wsize = 101;
opts.eltoldp = false;
opts.plot = false;
opts = parsepropval(opts,varargin{:});

% Nr. of samples in interval
K = opts.nrint;

% Nr. of parameters
np = size(p_lim,1);

% Convolutional kernel size
% ha wsize = 2*k --> wsize = 2*k+1
% ha wsize = 2*k+1 --> wsize = 2*k+1
ws = opts.wsize + mod(opts.wsize+1,2);

% Kernel window radius
wr = (ws-1)/2;

% Legyen K-val oszthato
wr = wr + K - mod(wr,K);


% fig1 = figure;
% fig2 = figure;

%%

N_cell = cell(np,1);
T_cell = cell(np,1);
t_cell = cell(np,1);
p_cell = cell(np,1);
dp_cell = cell(np,1);

for i = 1:np
    %%
    szorzo = 1;
    
    for proba = 1:opts.maxnrit
        %% [0] Meghatarozom, hogy milyen hosszu jelet generaljak
        % 
        % Ennyivel hosszabb idointervallumon generalom jelet. Igy a vegen ossze
        % lehet majd nyomni ugy, hogy az osszenyomott jel tspan-je [t0 t1] legyen
        % (legalabb). Az utolso tag: (opts.wsize - 1) / 10 egy heurisztika. A
        % simitastol is fugg ugyanis, hogy mennyire kicsi lesz a derivalt jel
        % maximuma.

        Time_multiplier = szorzo * min(abs(dp_lim(i,:))) * (opts.wsize - 1) / 5;
        T_ = T * Time_multiplier;
        
        Ni = round(T_ / dt / K);
        N = Ni * K;

        % Generating convolutional kernel (szandekosan felakkora, mivel
        % ketszer konvolvalok).
        w = hamming(wr+1);
        w = w / sum(w);

        %% [1] RANDOM GENERATOR

        % [1.1] Random Gaussian noise
        rgn = randn(1,N + 2*wr);
        rgn = rgn - mean(rgn);
        
        % [1.2] Eloszor dp-t generaljuk le (es simitjuk)
        dp = conv(rgn, w, 'valid');

        assert(numel(dp) == N + wr, 'numel(dp) = %d != %d = N+wr', numel(dp), N+wr)
        
        % [1.3.1] Ezzel azt garantalom, hogy p(i Ni) = 0, ahol i = 1..K
        dp = dp - reshape(ones(Ni+wr/K,1)*mean(reshape(dp,[Ni+wr/K,K]),1),[1,N+wr]);
        dp = conv(dp,w,'valid');
        
        assert(numel(dp) == N, 'numel(dp) = %d != %d = N', numel(dp), N)
        
        % [1.3.2] Ezzel azt garantalom, hogy p(1) = p(end) = 0
        % dp = dp - mean(dp);
        
        %% [2] Eltoljuk y iranyba dp-t, hogy min/max arany jo legyen
        % 
        % Ez lenyegeben az [1.3] szamolast felulirja, kiveve, ha:
        % dp_lim=[-a,a].
        %
        % Adott a,b, u.h. a+b = adott
        % 
        % Keresem c-t u.h. 
        % 
        %   a'     a - c
        %  ---- = ------- = p  --->  a - c = p b - p c
        %   b'     b - c           (p-1) c = p b - a
        % 
        % Tehat:
        % 
        %      a - p b
        %  c = -------
        %       p - 1
        % 

        if opts.eltoldp
            % Jelenlegi min/max (dp):
            a = min(dp);
            b = max(dp);

            % Elvart min/max arany:
            p = dp_lim(i,1)/dp_lim(i,2);

            % Eltolas mertekenek kiszamitasa a fenti keplet alapjan
            c = ( a - p*b ) / ( p - 1 );

            % Itt tortenik az y-eltolas
            dp = dp + c;
        end

        % Kiszamoljuk p-t (integralunk)
        p = cumsum(dp) * dt;

        %% [3] Skalazzuk p-t, hogy kitoltse p_lim-et

        % Ehhez kell min/max:
        minp = min(p);
        maxp = max(p);

        % Skalazo tenyezo:
        scale = (p_lim(i,2) - p_lim(i,1))/(maxp - minp);

        % Skalazzuk es y iranyba eltoljuk p-t
        p = (p - minp) * scale + p_lim(i,1);

        % Skalazzuk a derivaltjat is
        dp = dp * scale;

        %% [4] Ugy hatarozzuk meg az idoskalat, hogy dp kitoltse dp_lim-et
        % 
        % Alapja:
        % 
        %  d/dt( f(at) ) = [ df/dt ]_{t = at} * a
        % 

        a = min([ dp_lim(i,1) / min(dp), dp_lim(i,2) / max(dp), T_/T ]);

        T_ = T_ / a;

        N = numel(p);

        t = linspace(0,T_,N);
        dp = dp * a;

        ROSSZASAG = T_*max(dp)/(T*dp_lim(i,2));
        if 0.5 < ROSSZASAG && ROSSZASAG < 2
            break
        end

        szorzo = 1/ROSSZASAG;

    end % for
    
    dp_comp = diff(p(:,[1,1:end]),1,2)/(T_/N);
        
    T_cell{i} = T_;
    N_cell{i} = N;
    t_cell{i} = t';
    p_cell{i} = p';
    dp_cell{i} = dp';
    
    %{

    figure(fig1)
    ax(1) = subplot(2*np,1,i); plot(t,p), title(sprintf('p%d',i)), axis tight
    ax(2) = subplot(2*np,1,np+i); plot(t,[ dp ; dp_comp ]), title 'dp, dpComp', axis tight   
    linkaxes(ax,'x')

    %}
end % for

[ ~ , index ] = min([ T_cell{:} ]);

t = t_cell{index};

%%

p_cell = cellfun(@(tt,pp) {interp1(tt,pp,t)}, t_cell, p_cell);
dp_cell = cellfun(@(tt,dpp) {interp1(tt,dpp,t)}, t_cell, dp_cell);

p = horzcat(p_cell{:});
dp = horzcat(dp_cell{:});

pcz_dispFunctionEnd(TMP_bzqpLWoBcwMvZxPDPVtb);


if opts.plot
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

%{

T_ = t(end);
N = numel(t);
dp_comp = diff(p(:,[1,1:end]),1,2)/(T_/N);
figure(fig2);
subplot(311), plot(t,p), title p
subplot(312), plot(t,dp), title dp
subplot(313), plot(t,dp_comp), title dpComp

%}

end % function




function test1
%%

p_lim = [
    0 1
    -2 0
    ];

dp_lim = [
    -10 10
    -1 1
    ];

T = 10;
dt = 0.1;

[t,p,dp] = pcz_pdp_lim_smooth_dp(T,dt,p_lim,dp_lim,'eltoldp',false,'plot',true);
size(t)
size(p)
size(dp)

end


function test2
%%

p_lim = [
    0 1
    ];

dp_lim = [
    -1 1
    ];

T = 10;
dt = 0.1;

[t,p,dp] = pcz_pdp_lim_smooth_dp(T,dt,p_lim,dp_lim,'eltoldp',false,'plot',true);
xlim([0 T])

end