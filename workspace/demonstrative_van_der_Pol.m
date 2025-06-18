%% 
%  File: Lorenz_attractor.m
%  Directory: 7_ftools/workspace/07_DOA_positive_systems
%  Author: Peter Polcz (ppolcz@gmail.com) 
% 
%  Created on 2021. March 02. (2020b)

%%
% Automatically generated stuff

G_reset(01)
%       │└─ verbosity (0:false, 1:true, 2:do not change)
%       └── scope depth (0:reset, 1:reset if not set)

He = he;
I = @(varargin) eye(varargin{:});
O = @(varargin) zeros(varargin{:});

P_init(12.02)

%% model

x_lim_theoretical = [
    -2 2
    -2 2
    ];
[x_lfr,x_lfr_cell] = pcz_generateLFRStateVector('x',x_lim_theoretical);
syms t x1 x2 real
x = [x1;x2];

mu = 1;

f = -[
    x2
    mu*(1 - x1^2)*x2 - x1
    ];

Clip = [-1 1];
x_lim = [
    Clip;
    Clip;
    ];

f_ode = matlabFunction(f,'vars',{t,x});

%%
%{
    
    f_ode = matlabFunction(f,'vars',{t,x});
    
    [a,b,c] = ndgrid(linspace(0,c1,7));
    IC = [a(:),b(:),c(:)];
    
    fig = figure(141);
    delete(fig.Children);
    ax = axes(fig);
    hold on, grid on, box on;
    
    for x0 = IC'
    
        x0 = x0 + 0.1*randn(3,1);
        x0(x0 < 0) = 0;
    
        [~,xx] = ode45(f_ode,[0,1],x0);
        plot3(xx(:,1),xx(:,2),xx(:,3),'Color',[0.5 0.5 1]*0.9);
        plot3(xx(end,1),xx(end,2),xx(end,3),'r.')
        drawnow
    end

%}
%%
%{
    
C = 8;
[a,b] = ndgrid( ...
    linspace(x_lim(1,1),x_lim(1,2),15), ...
    linspace(x_lim(2,1),x_lim(2,2),15));
IC = [a(:),b(:),C - a(:) - b(:)];

fig = figure(141);
delete(fig.Children);
ax = axes(fig);
hold on, grid on, box on;

xEnd = zeros(3,size(IC,1));

for i = 1:size(IC,1)

    x0 = IC(i,:)';
    x0(x0 < 0) = 0;

    [~,xx] = ode45(dfx_ode,[0,1],x0);
    plot3(xx(:,1),xx(:,2),xx(:,3),'Color',[0.5 0.5 1]*0.9);
    plot3(xx(end,1),xx(end,2),xx(end,3),'r.')
    drawnow
    xEnd(:,i) = xx(end,:)';
end

%}
%%

f_fh_cell = cellfun(@(f) {matlabFunction(f,'vars',x)}, num2cell(dfx));
f_lfr_cell = cellfun(@(f) {f(x_lfr_cell{:})},f_fh_cell);

%%

x1 = x_lfr_cell{1};
x2 = x_lfr_cell{2};

Pi_lfr = [
    x1
    x2
    x1^2
    x1*x2
    x2^2
    x1^3
    x1^2*x2
    x1*x2^2
    x2^3
    ];
nBlock = numel(psi(0));

Pi = plfr(Pi_lfr);
m = size(Pi_lfr,1); % nBlock*3 = m

dPi = Pi.diff(x_lfr_cell,f_lfr_cell);

Pid0 = plfr([ Pi_lfr ; dPi.lfrtbx_obj ]);
J = [ I(m) O(m) ];
Jd = [ O(m) I(m) ];

Pid1 = Pid0.generatePI;

[Sd,Pid,iSd,Kerd] = P_mingen_for_LFR(Pid1);

Ed = J*[ Pid0.A Pid0.B ]*Sd;
Ad = Jd*[ Pid0.A Pid0.B ]*Sd;

N = P_affine_annihilator_for_LFR(Pi,x_lfr);
Nd = P_affine_annihilator_for_LFR(Pid,x_lfr);

pcz_fhzero_report(dPi - Ad*Pid,x,1e-7, 'dot pi(x) = Ad pi(x)');

[s,m] = size(N);
[sd,md] = size(Nd);

%%

[X_v,~,X_fci,proj] = P_ndnorms_of_X(x_lim);

X_NrV = size(X_v,1);
[X_NrF,X_NrFc] = size(X_fci);

P = sdpvar(m);
L = sdpvar(m,s,'full');
Ld = sdpvar(md,sd,'full');

% CONS = [ 0.00001 <= c2 , c2 <= c3*0.0199999];
CONS = [];

for i = 1:X_NrV
    xi = X_v(i,:)';
    
    CONS = [ CONS,
        P + He( L*N(xi) ) - 1e-7*eye(m) >= 0
        He( Ed'*P*Ad + Ld*Nd(xi) ) <= 0
        ];
end

% Boundary conditions
tau = sdpvar(1,X_NrF);
Qf = @(alpha) P + blkdiag(-alpha,O(m-1));
for i = 1:X_NrF

    Verbosity = G_VERBOSE(0);
    [S_Fi,Pi_Fi,iS_Fi,Ker_Fi] = P_mingen_for_LFR(Pi,'proj',proj{i});
    N_Fi = P_affine_annihilator_for_LFR(Pi_Fi,x_lfr,'proj',proj{i});
    G_VERBOSE(Verbosity);

    L_b1 = sdpvar(N_Fi.nu,N_Fi.ny);
    L_b2 = sdpvar(N_Fi.nu,N_Fi.ny);

    for j = 1:X_NrFc
        x_num = X_v(X_fci(i,j),:)';
        CONS = [ CONS
            S_Fi'*Qf(1)*S_Fi + He( L_b1 * N_Fi(x_num) ) >= 0
            % S_Fi'*Qf(tau(i))*S_Fi + He( L_b2 * N_Fi(x_num) ) <= 0
            ];
    end
end

sol = optimize(CONS,[],sdpsettings('solver','mosek'))
% sol = optimize(CONS,sum(tau),sdpsettings('solver','mosek'))
pcz_feasible(sol,CONS)


%%%

L = double(L);
P = double(P);

Pi_sym = sym(Pi);
V_fh = matlabFunction(Pi_sym.' * P * Pi_sym,'Vars',x);

[a,b] = ndgrid( ...
    linspace(x_lim(1,1),x_lim(1,2),81), ...
    linspace(x_lim(2,1),x_lim(2,2),81));
V = V_fh(a,b);

%%

fig = figure(1);
delete(fig.Children)
hold on, grid on, box on;

Logger.latexify_axis(gca,14)
Logger.latexified_labels(gca,14,'$x_1$','$x_2$')

surf(a,b,V)
contour(a,b,V)

x_limT = x_lim';
axis(x_limT(:)');
