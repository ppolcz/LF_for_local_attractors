%% local_brm_simplex_LFR
%
%  File: local_brm_simplex_LFR.m
%  Directory: 7_ftools/workspace/07_DOA_positive_systems
%  Author: Peter Polcz (ppolcz@gmail.com)
%
%  Created on 2019. December 30. (2019b)
%  Revised on 2024. October 10. (2024b)
%
% EXECUTIONS:
%  - 2021.03.02. (március  2, kedd), 12:49
%  - 2024.10.10. (October 10, Thursday), 13:56
% 
% │   ┌ local_brm_simplex_LFR:229 (pcz_fhzero_report) - projection on facet 4
% │   │   [  INFO  ] Tolerance: 1e-10.
% │   │   - Maximal difference: 2.98586e-11
% │   │   ----------------------------------------------------
% │   │   [   OK   ] projection on facet 4
% │   └ 0.006359 [sec]
% │    
% │   ┌ local_brm_simplex_LFR:232 (pcz_fhzero_report) - annihilator on facet 4
% │   │   [  INFO  ] Tolerance: 1e-10.
% │   │   - Maximal difference: 1.49707e-15
% │   │   ----------------------------------------------------
% │   │   [   OK   ] annihilator on facet 4
% │   └ 0.006033 [sec]
% │    
% │   [ FAILED ] Numerical problems (MOSEK). Solver time: 0.192755
% │   [ FAILED ] The solution is NOT feasible. Tolerance: 1e-10. Min(Primal) = -1.677484e+01.
% │    
% │   ┌ pcz_feasible_Finsler:18 - Positivity of P + He{L N}
% │   │   [  INFO  ] Mode: vertices are given, nr. of corners: 4
% │   │   [  INFO  ] Tolerance: 1e-10, positive tolerance: 1e-06.
% │   │   [  INFO  ] LMI size: (7x7), annihilator rows: 8. 
% │   │   ----------------------------------------------------
% │   │   [   OK   ] This LMI is feasible along the given tolerance value.
% │   └ 0.00953 [sec]
% │    
% │   ┌ pcz_feasible_Finsler:18 - Negativity of Pd + He{LdN d}
% │   │   [  INFO  ] Mode: vertices are given, nr. of corners: 4
% │   │   [  INFO  ] Tolerance: 1e-10, positive tolerance: 1e-06.
% │   │   [  INFO  ] LMI size: (17x17), annihilator rows: 22. 
% │   │   ----------------------------------------------------
% │   │   [   OK   ] This LMI is feasible along the given tolerance value.
% │   └ 0.004778 [sec]
% └ 15.4672 [sec]

%%
% Automatically generated stuff

G_reset(01)

try c = evalin('caller','persist'); catch; c = []; end
persist = Persist(mfilename('fullpath'), c); clear c;
persist.backup();
%clear persist

P_init(12.02)

He = he;
I = @(varargin) eye(varargin{:});
O = @(varargin) zeros(varargin{:});

%%

% 2021.03.02. (március  2, kedd), 12:45
x_lim_theoretical = [
    0 51
    0 51
    ];
[x_lfr,x_lfr_cell] = pcz_generateLFRStateVector('x',x_lim_theoretical);
pcz_generateSymStateVector(2,'x');
syms t real

V = 4;
Sf = 10;
Y = 0.5;
mu_max = 1;
K1 = 0.03;
K2 = 0.5;

%                       S
% μ(S) = μ_max * -----------------
%                 K2*S^2 + S + K1
%
mu = @(S) mu_max * S ./ (K2*S.^2 + S + K1);

% equilibrium point (S0,X0,F0)
S0 = 0.5 * (-2*K1 + 2*sqrt(K1^2+Sf^2*K1*K2 + Sf*K1)) / (Sf*K2 + 1);
X0 = (Sf - S0) * Y;
F0 = mu(S0)*V;

X = x1;
S = x2;

k = 1;

F = F0 - k*(S - S0);

f = [
    mu(S)*X - X*F/V
    -mu(S)*X/Y + (Sf-S)*F/V
    ];

f_fh = matlabFunction(f);

xeq = [ X0 ; S0 ];

%% Simulate on grid
%{

fig = figure(5);
ax = axes('Parent',fig); hold on

[xx1,xx2] = meshgrid(linspace(X0-1,X0+1,40),linspace(S0-0.1,S0+0.1,40));

xx1 = xx1(:);
xx2 = xx2(:);

f_ode = matlabFunction(f,'vars',{t x});

for i = 1:numel(xx1)
    [tt,xx] = ode45(f_ode,[0,4],[xx1(i) ; xx2(i)]);
    plot(xx(:,1),xx(:,2));
    plot(xx(1,1),xx(1,2),'.');
    drawnow
end

%}

%% Simulate in specific points
%{

fig = figure(5);
ax = axes('Parent',fig); hold on

xx0 = [
   0.0001 0.0001
   ];


f_ode = matlabFunction(f,'vars',{t x});

for x0 = xx0'
    x0
    [tt,xx] = ode45(f_ode,[0,100],x0);
    plot(xx(:,1),xx(:,2));
    plot(xx(1,1),xx(1,2),'.');
    drawnow
end

%}

%%

fi_fh = cellfun(@(fi) {matlabFunction(fi,'vars',x,'Optimize',false)},num2cell(f));
fi_lfr = cellfun(@(fi) {fi(x_lfr_cell{:})}, fi_fh);

%%

f_lfr = f_fh(x_lfr_cell{:});
f_plfr = plfr(f_lfr);

Pi_init = f_plfr.generatePI;

[S,Pi,iS,Ker] = P_mingen_for_LFR(Pi_init);

m = Pi.ny;
M = [ f_plfr.A f_plfr.B ] * S;

dPi = plfr(diff(Pi,x_lfr_cell,fi_lfr));

Pid0 = plfr([ Pi.lfrtbx_obj ; dPi.lfrtbx_obj ]);
J = [ I(m) O(m) ];
Jd = [ O(m) I(m) ];

Pid1 = Pid0.generatePI;

[Sd,Pid,iSd,Kerd] = P_mingen_for_LFR(Pid1);

Ed = J*[ Pid0.A Pid0.B ]*Sd;
Ad = Jd*[ Pid0.A Pid0.B ]*Sd;

% Ellenorzes
pcz_symzero_report(M*sym(Pi) - f, 'f(x) = A pi(x)')
% pcz_symzero_report(sym(dPi) - jacobian(sym(Pi),x)*f, 'differentiation with LFR vs. symbolic operations')
pcz_fhzero_report(plfr(Ed*Pid.lfrtbx_obj - Pi.lfrtbx_obj),x,1e-9, 'pi(x) = Ed pid(x)')
pcz_fhzero_report(plfr(Ad*Pid - dPi.lfrtbx_obj),x,1e-9, 'jacobian(pi(x)) f(x) = Ad pid(x)')

% Annihilators
N = P_affine_annihilator_for_LFR(Pi,x_lfr);
Nd = P_affine_annihilator_for_LFR(Pid,x_lfr);

pcz_fhzero_report(plfr(N * Pi),x,1e-10, 'N(x) pi(x) = 0')
pcz_fhzero_report(plfr(Nd * Pid),x, 1e-9, 'Nd(x) pid(x) = 0')

% Sizes
m = Pi.ny;
s = N.ny;
md = Pid.ny;
sd = Nd.ny;

%%

% Bounds
% x_val = [10 4];
% x_offset = [0 0];
% [X_v,X_fci,proj,~,in_simplex] = P_ndnorms_of_simplex(x_val,x_offset);

x_lim = [
    0 7
    0 3.5
    ];
[X_v,~,X_fci,proj] = P_ndnorms_of_X(x_lim);

X_NrV = size(X_v,1);
[X_NrF,X_NrFc] = size(X_fci);

P = sdpvar(m);
L = sdpvar(m,s,'full');
Ld = sdpvar(md,sd,'full');

CONS = [];

for i = 1:X_NrV
    xi = X_v(i,:)';

    CONS = [ CONS,
        P + He( L*N(xi) ) - 1e-5*I(m) >= 0
        He( Ed'*P*Ad + Ld*Nd(xi) ) + 0*I(md) <= 0
        ];
end

% Central condition
PiEq = Pi(xeq);
CONS = [ CONS , PiEq' * P * PiEq <= 0.1 ];

% Boundary conditions
tau = sdpvar(1,X_NrF);
Qf = @(alpha) P + blkdiag(-alpha,O(m-1));
for i = 1:X_NrF
%%
    [S_Fi,Pi_Fi,iS_Fi,Ker_Fi] = P_mingen_for_LFR(Pi,'proj',proj{i});
    N_Fi = P_affine_annihilator_for_LFR(Pi_Fi,x_lfr,'proj',proj{i});

    pcz_fhzero_report(@(x) Pi(proj{i}(x)) - S_Fi*Pi_Fi(proj{i}(x)),x,1e-10,...
        'projection on facet %d', i);

    pcz_fhzero_report(@(x) N_Fi(proj{i}(x))*Pi_Fi(proj{i}(x)),x,1e-10,...
        'annihilator on facet %d', i);

    L_b1 = sdpvar(N_Fi.nu,N_Fi.ny);
    L_b2 = sdpvar(N_Fi.nu,N_Fi.ny);

    for j = 1:X_NrFc
        x_num = X_v(X_fci(i,j),:)';
        CONS = [ CONS
            S_Fi'*Qf(1)*S_Fi + He( L_b1 * N_Fi(x_num) ) >= 0
            S_Fi'*Qf(tau(i))*S_Fi + He( L_b2 * N_Fi(x_num) ) <= 0
            ];
    end
end

OBJ = sum(tau);

sol = optimize(CONS,OBJ)
pcz_feasible(sol,CONS)

pcz_feasible_Finsler(P,L,N,X_v,'title','Positivity of P + He{L N}')
pcz_feasible_Finsler(-Ed'*P*Ad,-Ld,Nd,X_v,'title','Negativity of Pd + He{LdN d}')





%% Visualize

L = double(L);
Ld = double(Ld);
P = double(P);
Pi_sym = sym(Pi);

V = Pi_sym'*P*Pi_sym;
dV = matlabFunction(jacobian(V)*f);
V = matlabFunction(V);

% [xx1_HD,xx2_HD] = meshgrid(linspace(x_offset(1),x_val(1),61),linspace(x_offset(2),x_val(2),61));
[xx1_HD,xx2_HD] = meshgrid(linspace(x_lim(1,1),x_lim(1,2),161),logspace(-5,log10(x_lim(2,2)),161));
V_num_HD = V(xx1_HD,xx2_HD);

[xx1,xx2] = meshgrid(linspace(x_lim(1,1),x_lim(1,2),31),logspace(-2,log10(x_lim(2,2)+0.01),31)-0.01);
V_num = V(xx1,xx2);
dV_num = dV(xx1,xx2);

W = [ xx1(:) , xx2(:) ]';

% Compute contour
fig = figure;
C1 = pcontour(xx1_HD,xx2_HD,V_num_HD,[1 1]*1);
[~,C1_ind] = max([C1.nr]);
C1_cell = num2cell(C1(C1_ind).contour,2);
C1_cell{3} = C1_cell{1}*0 + min(V_num(:)) + eps;
close(fig);

threshold_V = 2.1;
threshold_dV = -20;
V_num(V_num > threshold_V) = nan;
dV_num(dV_num < threshold_dV) = nan;

X_ind = convhull(X_v(:,1),X_v(:,2));
X_v1 = X_v(X_ind,:);

fig = figure(3);

% Lyapunov function and contour
ax1 = subplot(121); hold off
T = delaunay(xx1,xx2);
TS_V = trisurf(T,xx1,xx2,V_num);
TS_V.LineWidth = 0.1;
TS_V.EdgeAlpha = 1;
TS_V.FaceColor = 'interp';


% shading interp,
% Lobj1 = light('Position',[-1 -1 4]);
% Lobj2 = light('Position',[-1 0.3 0]);
hold on
C1 = plot(C1_cell{1:2}, 'LineWidth',2);
X1 = plot(X_v1(:,1),X_v1(:,2),'o-');
Peq = plot(xeq(1),xeq(2),'.k','MarkerSize',20);
ax1.ZLim = [0 1.5];

% TS_dV2 = trisurf(T,xx1,xx2,dV_num);
% TS_dV2.CData = -V_num;

persist.latexify_axis(15)

Leg = plegend('$V(x)$','$\Omega$','$\mathcal X$','$x^*$');
Leg.Location = 'northwest';
Leg.FontSize = 15;

persist.latexified_labels(gca,18,'$x_1$','$x_2$');
ax1.XTick = 0:ceil(x_lim(1,2));
ax1.YTick = 0:ceil(x_lim(2,2));

% C1 = plot3(C1_cell{:}, 'LineWidth',2);
% X1 = plot3(X_v1(:,1),X_v1(:,2),X_v1(:,1)*0+min(V_num(:))+eps,'o-');

% Derivative
ax2 = subplot(122); hold off
TS_dV = trisurf(T,xx1,xx2,dV_num);
TS_dV.LineWidth = TS_V.LineWidth;
TS_dV.EdgeAlpha = TS_V.EdgeAlpha;
TS_dV.FaceColor = TS_V.FaceColor;
% ax2.ZLim = [-15 0];
ax2.XTick = ax1.XTick;
ax2.YTick = ax1.YTick;
hold on

persist.latexify_axis(15)

Leg = plegend('$\dot V(x)$');
Leg.Location = 'northeast';
Leg.FontSize = 15;

persist.latexified_labels(gca,18,'$x_1$','$x_2$');

colormap jet

%{

pngname = persist.file_simple('fig','bioreactor_rect.png');
print(pngname,'-dpng','-r500')

%}

persist.stoplog;

AREA_Omega = polyarea(C1_cell{1:2})
AREA_X = polyarea(X_v1(:,1),X_v1(:,2))

%%

function SYMBOLIC_MODEL_TO_LATEX
%%

syms x1 x2 V S_F Y mu K1 K2 X0 S0 F0 k real
syms mu(x2)

X = x1;
S = x2;

mu_expr = mu_max * S ./ (K2*S.^2 + S + K1);


F = F0 - k*(S - S0);

f = [
    mu(S)*X - X*F/V
    -mu(S)*X/Y + (S_F-S)*F/V
    ]

latex(f)

end

