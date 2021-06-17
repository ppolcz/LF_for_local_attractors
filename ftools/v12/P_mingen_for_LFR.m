function [S,syslfr_min,iS,Ker] = P_mingen_for_LFR(syslfr,varargin)
%% P_mingen_for_LFR
%  
%  File: P_mingen_for_LFR.m
%  Directory: 1_PhD_projects/00_my_toolboxes/FinslerTools/v11
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. March 11.
%

%% Initialization

args.State = [];
args.Round = 10;
args.proj = [];
args.lims = [];
args = parsepropval(args, varargin{:});

TMP_EOLHIWBuIwgwCuuNBzen = pcz_dispFunctionName;

%%

PI = syslfr;
if isa(PI,'lfr')
    PI = plfr(syslfr);
end
Pi = PI;

nu = PI.nu;
m = PI.ny;

if ~isempty(args.State) && isa(args.State,'sym')
    error 'TODO: EZ MEG NEM JO'
    % Nonlinear case: f(x,p) = A(x,p)*x = [F11 F12] * [x;pi_1(x,p)]
    %  where pi(x,p) = [x;pi_1(x,p)] = [I;PI_1(x,p)] * x

    x = args.State; %#ok<UNRCH>
    one = ones(size(x));
    x_lfr = plfr(one*0,diag(one),one,diag(one*0),diag(x),[-one one]);
  
    PI = Pi;
    Pi = plfr(Pi * x_lfr);
end

if ~isempty(args.State) && isa(args.State,'lfr')
    % Nonlinear case: f(x,p) = A(x,p)*x = [F11 F12] * [x;pi_1(x,p)]
    %  where pi(x,p) = [x;pi_1(x,p)] = [I;PI_1(x,p)] * x

    x = args.State;
  
    PI = Pi;
    Pi = plfr(Pi * x);
end

%%

rcond_tol = 1e-7;

%%

% 2020.05.20. (május 20, szerda), 01:52 (v12)
Ker = P_constKer_for_LFR(Pi','proj',args.proj,'lims',args.lims)';

Theta = null(Ker);

%% Find V and W, such that V is invertible

% Tolerance to rref (for rank decision).
rank_tol = max(size(Theta)) * eps(1); % norm(Theta) = 1
while true

    [~,Irows] = rref(Theta', rank_tol);
    V = Theta(Irows,:);

    if rcond(V*V') < rcond_tol
        % sigma = svd(Theta)';
        % pcz_dispFunction_num2str(sigma)
        
        rank_tol = rank_tol * 10;
        pcz_dispFunction('Rank decision tolerance for `rref'' increased to %g', rank_tol)
        continue
    end
    
    break
end


Jrows = setdiff(1:m,Irows);
W = Theta(Jrows,:);

sigma = [Irows Jrows];
Is = pcz_permat(sigma)';

% pcz_display(Theta,V,W,Irows,Jrows,sigma)

%% Generate transformation matrix

Gamma = W*V'/(V*V');

%{
    args.Round = 10;
%}

if args.Round
    Gamma = round(Gamma, args.Round);
end

% Model reduction transformation driven by:
% I_sigma*PI = S_sigma*wh_PI
% Ss: S_sigma
Ss = [
    eye(size(Gamma,2))
    Gamma  
    ];

S = Is'*Ss;

iS = [ eye(size(Gamma,2)) 0*Gamma' ] * Is;

syslfr_min = plfr(iS * PI);
% syslfr_min = syslfr_min.set_vars(PI.subsvars)

%{    

    sym(Pi) - S*sym(syslfr_min)

    pcz_display(V,W,Gamma,S,pInv_S)
    pcz_display(Is_left,M,Is_right,Ss)

%}


pcz_dispFunctionEnd(TMP_EOLHIWBuIwgwCuuNBzen);

%% Igy lehetne sokkal rovidebben
% 
% Ker = P_constKer_for_LFR(syslfr')';
% 
% S = null(Ker);
% 
% iS = pinv(S);
% 
% syslfr_min = minlfr(iS*syslfr);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function test1_simple
%%
% 2020.09.01. (szeptember  1, kedd), 10:10

x = pcz_generateLFRStateVector('x',[-1 1 ; -2 3 ; -3 4],[-1 1 ; -2 2 ; -3 3]);

x1 = x(1);
x2 = x(2);
x3 = x(3);

Pi_init = [
    1
    x
    x1+x1*x2
    x1*x2
    x1^2*x2
    ];

[S,Pi,is,Ker] = P_mingen_for_LFR(Pi_init)

end

function test2_simple_proj
%%
% 2020.09.01. (szeptember  1, kedd), 10:28

x_lim = [
    -1 1 
    -2 3 
    -3 4
    ];

p_lim = [
    3 6
    0 1
    ];

dp_lim = [
    -1 1
    -2 2
    ];

x = pcz_generateLFRStateVector('x',x_lim);
x1 = x(1);
x2 = x(2);
x3 = x(3);

p = pcz_generateLFRStateVector('p',p_lim,dp_lim);
p1 = p(1);
p2 = p(2);

[X_v, X_fci, X_proj_x, X_proj_xp, X_in_simplex, X_Vk] = P_ndnorms_of_simplex([1 2 3],[0.1 0.1 0.1],p_lim)


X_proj_xp{4}([-1;2;4;4;1])


Pi_init = [
    1
    x
    p1
    x*p1
    x*p2
    x1+x1*x2
    x1*x2
    x1^2*x2
    ];

[S,Pi,is,Ker] = P_mingen_for_LFR(Pi_init,'proj',X_proj_xp{4},'lims',[x_lim ; p_lim])

sym(Pi)

end