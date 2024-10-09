%% 28. Genetic toggle switch [VEDETT]
%
%  Author: Peter Polcz (ppolcz@gmail.com)
%  Created on 2020. January 04. (2019b)
%  Minor revision on 2021. June 17. (2021a)
%
% EXECUTIONS:
% - 2024.10.03. (October  3, Thursday), 11:14 (Matlab R2023b)

%%
% Automatically generated stuff

G_reset(01)
%       |L- verbosity (0:false, 1:true, 2:do not change)
%       L-- scope depth (0:reset, 1:reset if not set)

He = he;
I = @(varargin) eye(varargin{:});
O = @(varargin) zeros(varargin{:});

% 2021.03.02. (mrcius  2, kedd), 09:42
P_init(12.01)

%%

% 2020.09.01. (szeptember  1, kedd), 11:19
x_lim_theoretical = [
    -10 10
    -10 10
    ];
[x_lfr,x_lfr_cell] = pcz_generateLFRStateVector('x',x_lim_theoretical);
pcz_generateSymStateVector(2,'x');
syms t real

alpha_1 = 1.3;
alpha_2 = 1;
beta = 3;
gamma = 10;

f = [
    alpha_1 / (1 + x2^beta) - x1
    alpha_2 / (1 + x1^gamma) - x2
    ];

f_ode = matlabFunction(f,'vars',{t, x});

xeq = [
    0.668   0.9829       % locally asymptotically stable
    0.8807  0.7808       % separatrix
    1.2996  0.0678       % locally asymptotically stable
    ];

%%

fi_fh = cellfun(@(fi) {matlabFunction(fi,'vars',x,'Optimize',false)},num2cell(f));
fi_lfr = cellfun(@(fi) {fi(x_lfr_cell{:})}, fi_fh);

%%

%{

f_ode = matlabFunction(f,'vars',{t x});
figure, hold on
for i = 1:100
    [tt,xx] = ode45(f_ode, [0,40],abs(20*randn(n_x,1)));
    P = plot(xx(:,1),xx(:,2),'.-')
    plot(xx(1,1),xx(1,2),'.','MarkerSize',10,'Color',P.Color)
end

%}

f_fh = matlabFunction(f,'vars',x,'Optimize',false);

f_lfr = f_fh(x_lfr_cell{:});
f_plfr = plfr(f_lfr);

%{
%% Ez mar konkretan a cikkhez kell

f11 = f_plfr.A;
F12 = f_plfr.B;
f21 = f_plfr.C;
F22 = f_plfr.D;
Delta = f_plfr.Delta;
pi_ = sym(f_plfr.generatePI);
m_ = size(F22,1);

pcz_num2str_latex(f11);
pcz_num2str_latex(F12);
pcz_num2str_latex(f21);
pcz_num2str_latex(F22);
pcz_latex_v7(Delta);
pcz_latex_v7(pi_);

%%
%}

Pi_init = f_plfr.generatePI;

[S,Pi,iS,Ker] = P_mingen_for_LFR(Pi_init);

%{
%% Ez mar konkretan a cikkhez kell

wh_pi_ = sym(Pi);
m_prime = size(F22,1);
pcz_num2str_latex(S);
pcz_latex_v7(wh_pi_);

%%
%}

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

%{
%% Ez mar konkretan a cikkhez kell

h11 = Pid0.A;
H12 = Pid0.B;
h21 = Pid0.C;
H22 = Pid0.D;
Delta_h = Pid0.Delta;

Sg = Sd;
% wh_pi_d = sym(Pid);
m_prime = size(F22,1);

pcz_num2str_latex(h11);
pcz_num2str_latex(H12);
pcz_num2str_latex(h21);
pcz_num2str_latex(H22);
pcz_latex_v7(Delta_h);
pcz_num2str_latex(Sg);
% pcz_latex_v7(wh_pi_);

g11 = Pid1.A;
G12 = Pid1.B;
g21 = Pid1.C;
G22 = Pid1.D;
Delta_g = Pid1.Delta; 
pcz_num2str_latex(g11);
pcz_num2str_latex(G12);
pcz_num2str_latex(g21);
pcz_num2str_latex(G22);
pcz_latex_v7(Delta_g);


A_ket_delta_ugyanaz = norm(double(Delta_g - Delta_h)) == 0

%%
%}



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

x_lim = [
    0.5 1.5
    -0.18 1.2
    ];

% Kiszamolom, hogy (O,E2) hol metszi X-et. Egyenes egyenlete:
E2_x = xeq(2,1);
E2_y = xeq(2,2);

% Baloldali metszespont
y_inters = x_lim(1,1) * E2_y / E2_x;
if ( x_lim(2,1) <= y_inters ) && ( y_inters <= x_lim(2,2) )
    P1 = [
        x_lim(1,1)
        y_inters
        ]';
else
    P1 = [
        x_lim(2,1) * E2_x / E2_y
        x_lim(2,1)
        ]';
end

% Jobboldali metszespont
y_inters = x_lim(1,2) * E2_y / E2_x;
if ( x_lim(2,1) <= y_inters ) && ( y_inters <= x_lim(2,2) )
    P2 = [
        x_lim(1,2)
        y_inters
        ]';
else
    P2 = [
        x_lim(2,2) * E2_x / E2_y
        x_lim(2,2)
        ]';
end

[X_v,~,X_fci,proj] = P_ndnorms_of_X(x_lim,'Vertices',[P1 ; P2],'Facets',[ 1 2 ]);
X_NrV = size(X_v,1);
[X_NrF,X_NrFc] = size(X_fci);


P = sdpvar(m);
L = sdpvar(m,s,'full');
Ld = sdpvar(md,sd,'full');
tau = sdpvar(1,X_NrF);

CONS = [
    Pi(xeq(1,:))' * P * Pi(xeq(1,:)) <= 0.3
    Pi(xeq(3,:))' * P * Pi(xeq(3,:)) <= 0.3
    ];

% OBJ = sum(tau);
OBJ = [];

for i = 1:X_NrV
    xi = X_v(i,:)';

    CONS = [ CONS,
        P + He( L*N(xi) ) - 1e-5*I(m) >= 0
        He( Ed'*P*Ad + Ld*Nd(xi) ) + 0*I(md) <= 0
        ];
end

% Boundary conditions
Qf = @(alpha) P + blkdiag(-alpha,O(m-1));
for i = 1:X_NrF

    % 2021.03.08. (mrcius  8, htf), 16:33
    % Szimbolikusan is mkdik, de gyorsabb numerikusan
    
    Verbosity = G_VERBOSE(0);
    [S_Fi,Pi_Fi,iS_Fi,Ker_Fi] = P_mingen_for_LFR(Pi,'proj',proj{i});
    N_Fi = P_affine_annihilator_for_LFR(Pi_Fi,x_lfr,'proj',proj{i});
    % N_Fi = P_affine_annihilator_for_LFR(Pi_Fi,x_lfr,'proj',proj{i},'sym',1,'lims',x_lim);
    G_VERBOSE(Verbosity);

    %{
    %% Ez mar konkretan a cikkhez kell

    pcz_num2str_latex(S_Fi);
    pcz_latex_v7(sym(Pi_Fi),'label',sprintf('\\\\pi_{f,%d} = ',i))
    pcz_latex_v7(sym(N_Fi),'label',sprintf('N_{f,%d} = ',i))

    %%
    %}
    
    L_b1 = sdpvar(N_Fi.nu,N_Fi.ny);
    L_b2 = sdpvar(N_Fi.nu,N_Fi.ny);

%{
    % Ellenorzeskeppen
    Ker_Fi = round(rref(Ker_Fi),10)
    N_Fi = P_affine_annihilator_for_LFR(Pi_Fi,x_lfr,'proj',proj{i},'sym',1);
    Pi_Fi_sym = sym(Pi_Fi)
    N_Fi_sym = sym(N_Fi)
    ZERO_elm = vpa(simplify(N_Fi_sym * Pi_Fi_sym),2)
%}

    for j = 1:X_NrFc

        x_num = X_v(X_fci(i,j),:)';

        CONS = [ CONS

            S_Fi'*Qf(1)*S_Fi + He( L_b1 * N_Fi(x_num) ) >= 0

            S_Fi'*Qf(tau(i))*S_Fi + He( L_b2 * N_Fi(x_num) ) <= 0

            ];

    end

end

sol = optimize(CONS,OBJ)
pcz_feasible(sol,CONS)

%% Visualize

L = double(L);
Ld = double(Ld);
P = double(P);
Pi_sym = sym(Pi);

% %{
%% Ez mar konkretan a cikkhez kell

pcz_num2str_latex(P);

%%
%}

[xx1,xx2] = meshgrid(linspace(x_lim(1,1),x_lim(1,2),41),linspace(x_lim(2,1),x_lim(2,2),41));

V = Pi_sym'*P*Pi_sym;
dV_fh = matlabFunction(jacobian(V)*f);
V_fh = matlabFunction(V);

V_num = V_fh(xx1,xx2);
dV_num = dV_fh(xx1,xx2);

% Generate grid
resolution = { 471 471 }.';
lscell = cellfun(@(c) {linspace(c{:})}, num2cell([ num2cell(x_lim) resolution ],2) );
xxfun = cell(1,3);
[xxfun{1:2}] = ndgrid(lscell{:});

% Retain invariant level set
figtmp = figure;
xxfun{3} = V_fh(xxfun{1:2});
Cont1 = pcontour(xxfun{:},[1 1]);
close(figtmp);

figtmp = figure;
xxfun{3} = V_fh(xxfun{1:2});
Cont2 = pcontour(xxfun{:},[0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9]);
close(figtmp);


%% Fancy plot a Lyapunov fuggveny

fig = figure(7);
delete(fig.Children);
ax = axes('Parent',fig);

Level_2D = -50;
Threshold = 100;

V_clip = V_num;
V_clip(V_clip > Threshold*1.4) = nan;

figure(fig);
hold off;
S2 = surf(xx1,xx2,V_clip);
S2.CData = log(V_clip);
hold on

x0s = [
    1.5 1.2
%     1.4 1.2
    1.3 1.2
%     1.2 1.2
    1.1 1.2
%     1.0 1.2
    0.9 1.2
%     0.8 1.2
    0.7 1.2
%     0.6 1.2
    0.5 1.2
%     1.5 -0.1
    1.3 -0.1
%     1.2 -0.1
    1.1 -0.1
%     1.0 -0.1
    0.9 -0.1
%     0.8 -0.1
    0.7 -0.1
%     0.6 -0.1
    0.5 -0.1
%     0.5 1.2
    0.5 1.1
%     0.5 1.0
    0.5 0.9
%     0.5 0.8
    0.5 0.7
%     0.5 0.6
    0.5 0.5
%     0.5 0.4
    0.5 0.3
%     0.5 0.2
    0.5 0.1
%     0.5 0.0
    0.5 -0.1
%     1.5 1.2
%     1.5 1.1
    1.5 1.0
%     1.5 0.9
    1.5 0.8
%     1.5 0.7
    1.5 0.6
%     1.5 0.5
    1.5 0.4
%     1.5 0.3
    1.5 0.2
%     1.5 0.1
    1.5 0.0
%     1.5 -0.1
    ]';
for x0 = x0s
    [tt,xx] = ode45(f_ode,[0 10],x0);
    % Pl1 = plot3(xx(:,1),xx(:,2),xx(:,2)*0+Level_2D,'Color',[1 1 1]*0.4);
end


for C = Cont1
    c = C.contour;
    plot3(c(1,:),c(2,:),c(2,:)*0+Level_2D,'LineWidth',2,'Color',pcz_get_plot_colors([],1))
end

% for i = 1:2
%     c = Cont(i).contour;
%     c = c(:,all(isfinite(c)))';
%     Fill = fill3(c(:,1),c(:,2),c(2,:)*0-2,pcz_get_plot_colors([],1)*1.3);
%     Fill.EdgeAlpha = 0;
%     Fill.FaceAlpha = 0.4;
% end

figure(fig);
for C = Cont2
    c = C.contour;
    plot3(c(1,:),c(2,:),c(2,:)*0+Level_2D,'LineWidth',0.5,'Color',pcz_get_plot_colors([],1))
    hold on
end

plot3(xeq(:,1),xeq(:,2),xeq(:,2)*0+Level_2D,'.k','MarkerSize',20)

axis(reshape(x_lim',[1 4]))
grid on

ax.View = [ -155.8045   27.7281 ];

ax.XTick = [x_lim(1,1):0.2:x_lim(1,2)];
ax.YTick = [x_lim(2,1) 0:0.2:1 x_lim(2,2)];
ax.ZTick = 0:20:Threshold;
ax.ZLim = [Level_2D Threshold];

axis vis3d
colormap jet

Logger.latexify_axis(12);
Logger.latexified_labels(ax,18,'$x_1$','$x_2$');

Mid1 = plot3([P1(1) P2(1)],[P1(2) P2(2)],[0 0]+Level_2D,'o--','Color',pcz_get_plot_colors([],5),'LineWidth',1.5);

X_ind = convhull(X_v(1:4,1),X_v(1:4,2));
X_v1 = X_v(X_ind,:);
X1 = plot3(X_v1(:,1),X_v1(:,2),X_v1(:,1)*0+Level_2D,'o-','Color',pcz_get_plot_colors([],2),'LineWidth',1.5);


Text1 = text(1.1,1,Level_2D,'$\mathcal{F}_5$',PersistStandalone.font_latex14{:});
Text1.VerticalAlignment = 'bottom';
Text1.HorizontalAlignment = 'right';

Text2 = text(1.45,0.6,Level_2D,'$\mathcal{X}$',PersistStandalone.font_latex14{:});

Text3 = text(0.7,0.8,Level_2D,'$x_*^{(1)}$',PersistStandalone.font_latex14{:});

Text4 = text(0.9,0.6,Level_2D,'$x_*^{(2)}$',PersistStandalone.font_latex14{:});

Text5 = text(1.3,0.4,Level_2D,'$x_*^{(3)}$',PersistStandalone.font_latex14{:});

Dirname = '/home/ppolcz/T/Dropbox/Peti/Munka/01_PPKE_2020/Dokumentaciok/28_nonzero_attractor_actual/fig/Exported_from_Matlab';
Figname = [ Dirname filesep datestr(today,29) '-genetic_toggle_switch' ];

%{

pause(1)
print([Figname '.png'],'-dpng','-r900')

pngname = persist.file_simple('fig','genetic_toggle_switch_intcurves.png');
print(pngname,'-dpng','-r900')

%}


