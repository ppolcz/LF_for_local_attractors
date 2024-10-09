%% local_Motzkin_system [VEDETT]
%
%  Author: Peter Polcz (ppolcz@gmail.com)
%  Created on 2019. November 11. (2019a)
%
% EXECUTIONS:
% - 2024.10.03. (October  3, Thursday), 11:14 (Matlab R2023b)

P_init(12)

%%

fname = pcz_mfilename;
cd(fname.dir)

savename = sprintf('../results/local_Motzkin_system_data-%s.mat',datestr(date,29));
save(savename,'fname')

%%

G_reset(01)
%        verbosity (0:false, 1:true, 2:do not change)
%        scope depth (0:reset, 1:reset if not set)

He = he;
I = @(varargin) eye(varargin{:});
O = @(varargin) zeros(varargin{:});

%%

% 2021.03.02. (mrcius  2, kedd), 12:15
x_lim = [
    -1 1
    -1 1
    ] * 2;
[x_lfr,x_lfr_cell] = pcz_generateLFRStateVector('x',x_lim);
pcz_generateSymStateVector(2,'x');
syms t real

H = x1^4*x2^2 + x1^2*x2^4 - 3*x1^2*x2^2 + 1;

% ezsurf(H,[-1.3,1.3],[-1.3,1.3])

f = -gradient(H,x);

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

% %{
%% Ez mar konkretan a cikkhez kell

f11 = f_plfr.A;
F12 = f_plfr.B;
f21 = f_plfr.C;
F22 = f_plfr.D;
Delta = f_plfr.Delta;
pi_sym = sym(f_plfr.generatePI);
m_ = size(F22,1);

pcz_num2str_latex(f11);
pcz_num2str_latex(F12);
pcz_num2str_latex(f21);
pcz_num2str_latex(F22);
pcz_latex_v7(Delta);
pcz_latex_v7(pi_sym);

f_lfr = f_plfr;
save_variables(savename,...
    f_lfr,'Dynamical equations $f(x)$',...
    pi_sym,'$\pi(x)$ originating from the dynamic equations, namely: $f(x) = A \pi(x)$')

%%
%}

Pi_init = f_plfr.generatePI;

[S,Pi,iS,Ker] = P_mingen_for_LFR(Pi_init);

% %{
%% Ez mar konkretan a cikkhez kell

pi_min_lfr = Pi;
pi_min_sym = sym(Pi);
m_prime = size(F22,1);
pcz_num2str_latex(S);
pcz_latex_v7(pi_min_sym);

save_variables(savename,...
    pi_min_lfr,'Dynamical equations $f(x)$',...
    pi_min_sym,'$\pi(x)$ originating from the dynamic equations, namely: $f(x) = A \pi(x)$')


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

% %{
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

% %{
%% Ez mar konkretan a cikkhez kell

N_sym = sym(P_affine_annihilator_for_LFR(Pi,x_lfr,'sym',true));
pcz_latex_v7(N_sym)


% Nd_affmat = paffmat(Nd,x);
% pcz_num2str_latex_output(Nd_affmat.Theta')

%%
%}

% 2019.12.28. (december 28, szombat), 02:40
[Y_v,~,Y_fci,Y_proj] = P_ndnorms_of_X(x_lim * blkdiag(0,1));
Y_v = [ Y_v(2:3,:) ; -Y_v(2:3,:) ];
Y_fci = [
    1 3
    2 4
    ];
Y_proj = Y_proj([1 3]);


[X_v,~,X_fci,proj] = P_ndnorms_of_X(x_lim);
X_NrV = size(X_v,1);

% 2019.12.28. (december 28, szombat), 02:41
X_v = [
    X_v
    Y_v
    ];
X_fci = [
    X_fci
    Y_fci + X_NrV
    ];
proj = [
    proj
    Y_proj
    ];
[X_NrF,X_NrFc] = size(X_fci);


P = sdpvar(m);
L = sdpvar(m,s,'full');
Ld = sdpvar(md,sd,'full');
tau = sdpvar(1,X_NrF);

% PiEq = Pi(xeq);

CONS = [
    Pi(+1,+1)' * P * Pi(+1,+1) <= 0.3
    Pi(+1,-1)' * P * Pi(+1,-1) <= 0.3
    Pi(-1,+1)' * P * Pi(-1,+1) <= 0.3
    Pi(-1,-1)' * P * Pi(-1,-1) <= 0.3
    ];
% CONS = [ ];

OBJ = sum(tau);
% OBJ = [];

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

    Verbosity = G_VERBOSE(0);
    [S_Fi,Pi_Fi,iS_Fi,Ker_Fi] = P_mingen_for_LFR(Pi,'proj',proj{i});
    N_Fi = P_affine_annihilator_for_LFR(Pi_Fi,x_lfr,'proj',proj{i},'sym',1,'lims',x_lim);
    G_VERBOSE(Verbosity);

    % %{
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

[xx1,xx2] = meshgrid(linspace(x_lim(1,1),x_lim(1,2),61),linspace(x_lim(1,1),x_lim(1,2),61));

V = Pi_sym'*P*Pi_sym;
dV_fh = matlabFunction(jacobian(V)*f);
V_fh = matlabFunction(V);

V_num = V_fh(xx1,xx2);
dV_num = dV_fh(xx1,xx2);

% Compute contour
fig = figure(5);

% Lyapunov function and contour
% S1 = surf(xx1,xx2,V_num); hold on
% C1 = plot3(C1_cell{:}, 'LineWidth',2);
contour(xx1,xx2,V_num,linspace(0,2,50));

V_num(V_num > 10) = nan;
dV_num(dV_num < -1) = nan;

% Lyapunov function
ax3 = figure(2);
S2 = surf(xx1,xx2,V_num);

% Derivative
ax3 = figure(3);
S2 = surf(xx1,xx2,dV_num);

%% Fancy plot a Lyapunov fuggveny vs Motzkin polinom

fig = figure(6);

H_fh = matlabFunction(H);
fh = { H_fh V_fh };
name = { '$H(x_1,x_2)$' '$V(x_1,x_2)$' };
alpha = { 0.99 1 };

shading_flat_or_light = 1;

if shading_flat_or_light
    resolution = { 51 51 }.';
else
    resolution = { 171 171 }.';
end

lscell = cellfun(@(c) {linspace(c{:})}, num2cell([ num2cell([-1 1 ; -1 1]*2) resolution ],2) );
xxfun = cell(1,3);
[xxfun{1:2}] = ndgrid(lscell{:});

threshold = 2.2;

for i = 1:2

    % Retain invariant level set
    figtmp = figure;
    xxfun{3} = fh{i}(xxfun{1:2});
    Cont = pcontour(xxfun{:},[1 1]*alpha{i});
    close(figtmp);

    figure(fig);
    ax = subplot(1,2,i); hold off;

    for C = Cont
        c = C.contour;
        plot3(c(1,:),c(2,:),c(2,:)*0,'LineWidth',2,'Color',pcz_get_plot_colors([],1))
        hold on
    end
    
    c = Cont(1).contour;
    c = c(:,all(isfinite(c)))';
    Fill = fill(c(:,1),c(:,2),pcz_get_plot_colors([],1)*1.3);
    Fill.EdgeAlpha = 0;
    Fill.FaceAlpha = 0.4;
    
    xxfun_mesh = cellfun(@(v) {v(xxfun{3} <= threshold)}, xxfun);
    tri = delaunay(xxfun_mesh{1:2});

    TS = trisurf(tri,xxfun_mesh{:});
    TS.LineWidth = 0.1;
    TS.FaceColor = 'interp';
    
    grid on
    ax.View = [  -71.1029   22.2111  ];
    ax.ZLim = [0 2];
    
    if ~shading_flat_or_light
        TS.EdgeAlpha = 0;
    
        Lgh1 = light;
        Lgh1.Position = [ -1 0 0 ];

        % Lgh2 = light;
        % Lgh2.Position = [ -1 1 0.3 ];
    end

    Logger.latexify_axis(15);
    Logger.latexified_labels(ax,22,'$x_1$','$x_2$');

    plot([-1 1 -1 1],[-1 -1 1 1],'.k','MarkerSize',20)
    
    colormap jet

    AREA_Gamma = polyarea(c(:,1),c(:,2));
    
    display(AREA_Gamma,sprintf('Area of sub-level set of %s', name{i}))
    
end


X_ind = convhull(X_v(1:4,1),X_v(1:4,2));
X_v1 = X_v(X_ind,:);
X1 = plot(X_v1(:,1),X_v1(:,2),'o-','Color',pcz_get_plot_colors([],2),'LineWidth',1.5);

MidX1 = plot([2 -2],[0 0],'o--','Color',pcz_get_plot_colors([],5),'LineWidth',1.5);
MidX2 = plot([0 0],[2 -2],'o--','Color',pcz_get_plot_colors([],5),'LineWidth',1.5);

Text1 = text(-2.25,-0.55,0.2,'$\Omega^{(1)}$',PersistStandalone.font_latex14{:});
Text2 = text(-1.6,-1.35,0.2,'$x_*^{(1)}$',PersistStandalone.font_latex14{:});

%{

pngname = persist.file_simple('fig','Motzkin_vs_LF_v2_flat.png');
print(pngname,'-dpng','-r500')

%}

%% Fancy plot a Motzkin polinomrol

fig = figure(6);

resolution = { 81 81 }.';
lscell = cellfun(@(c) {linspace(c{:})}, num2cell([ num2cell([-1 1 ; -1 1]*2) resolution ],2) );
xxfun = cell(1,3);
[xxfun{1:2}] = ndgrid(lscell{:});

H_fh = matlabFunction(H);
xxfun{3} = H_fh(xxfun{1:2});

xxfun = cellfun(@(v) {v(xxfun{3} <= 2)}, xxfun);

tri = delaunay(xxfun{1:2});

TS = trisurf(tri,xxfun{:});
TS.LineWidth = 0.1;
TS.EdgeAlpha = 0.3;
TS.FaceColor = 'interp';

ax = gca;

ax.View = [  -71.1029   22.2111  ];

Logger.latexify_axis(15);
Logger.latexified_labels(ax,22,'$x_1$','$x_2$','$H(x_1,x_2)$');

colormap jet

%{

pngname = persist.file_simple('fig','Motzkin_polynomial.png');
print(pngname,'-dpng','-r500')

%}

return 

%% Ezt sietve csak a command-line-ba rtam, szedd ssze ami rdekes:

contour(xx1,xx2,V_num,[1 1])
contour(xx1,xx2,V_num,[1 2])
contour(xx1,xx2,V_num,1:0.1:2)
min(V_num(:))
V_num - min(V_num(:))
xx1(31,31)
xx1(31,31),xx2(31,31)
V0_num = V_num - min(V_num(:))
V0_num(31,31)
V01_num = V0_num/V0_num(31,31)
V01_num
surf(xx1,xx2,V01_num)
surf(xx1,xx2,log(V01_num+1))
H
lH = matlabFunction(log(H+1))
figure
surf(xx1,xx2,lH(xx1,xx2))
figure
surf(xx1,xx2,lH(xx1,xx2)-log(V01_num+1))
H_fh = matlabFunction(H)
H_num = H_fh(xx1,xx2)
H_num(31,31)
min(H_num(:))
min(V_num(:))
V0_num = V_num - min(V_num(:)) + min(H_num(:))
min(V0_num(:))
V0_num(31,31)
V01_num = V0_num/V0_num(31,31)
V01_num(31,31)
min(V01_num(:))
V01_num = V01_num - min(V01_num(:)) + min(H_num(:))
V01_num = V01_num/V01_num(31,31)
V01_num = V01_num - min(V01_num(:)) + min(H_num(:))
V01_num = V01_num/V01_num(31,31)
V01_num = V01_num - min(V01_num(:)) + min(H_num(:))
V01_num = V01_num/V01_num(31,31)
V01_num = V01_num - min(V01_num(:)) + min(H_num(:))
V01_num = V01_num/V01_num(31,31)
min(V0_num(:))
V01_num(31,31)
surf(xx1,xx2,V01_num), hold on, surf(xx1,xx2,H_num)
surf(xx1,xx2,V01_num), hold on, surf(xx1,xx2,H_num), zlim([0 10])
V01n = V01_num;
V01n(V01n > 10) = nan
surf(xx1,xx2,V01n)
Hn = H_num
Hn(Hn > 10) = nan
hold on, surf(xx1,xx2,Hn)
contour(xx1,xx2,V_num,1:0.1:2)
contour(xx1,xx2,V_num,[1 1])
hold on, plot([-1 1 1 -1],[1 1 -1 -1],'k.','MarkerSize',10)
hold on, plot([-1 1 1 -1],[1 1 -1 -1],'k.','MarkerSize',30)
xlim(-2,2)
xlim([-2,2])
ylim([-2,2])
axis(equal)
axis('equal')
xlim([-2,2])
ylim([-2,2])
figure, subplot(121), surf(xx1,xx2,log(V01_num+1)), title('log(V(x) + 1)'), subplot(122),
surf(xx1,xx2,lH(xx1,xx2)), title('log(H(x) + 1)')
figure, ax(1) = subplot(121), surf(xx1,xx2,log(V01_num+1)), title('log(V(x) + 1)'), ax(2) = subplot(122), surf(xx1,xx2,lH(xx1,xx2)), title('log(H(x) + 1)'), linkaxes(ax)
figure, ax(1) = subplot(121), surf(xx1,xx2,log(V01_num+1)), title('log(V(x) + 1)'), ax(2) = subplot(122), surf(xx1,xx2,lH(xx1,xx2)), title('log(H(x) + 1)'), linkprop(ax, {'View', 'XLim', 'YLim', 'ZLim',})
shading interp
axis tight
figure, contour(xx1,xx2,H_fh(xx1,xx2),[0.9 1 1.1])
figure, contour(xx1,xx2,H_fh(xx1,xx2),[1 1.1])
figure, contour(xx1,xx2,H_fh(xx1,xx2),[0.99])
figure, contour(xx1,xx2,H_fh(xx1,xx2),[0.99 0.99])
hold on, plot([-1 1 1 -1],[1 1 -1 -1],'k.','MarkerSize',30)
figure, contour(xx1,xx2,H_fh(xx1,xx2),[0.6 0.7 0.8 0.9 0.99])
hold on, plot([-1 1 1 -1],[1 1 -1 -1],'k.','MarkerSize',30)
axis equal
xlim([-2,2])
ylim([-2,2])


figure, contour(xx1,xx2,V_num,[0.35 1 2])
hold on, plot([-1 1 1 -1],[1 1 -1 -1],'k.','MarkerSize',10)
hold on, plot([-1 1 1 -1],[1 1 -1 -1],'k.','MarkerSize',30)
xlim([-2,2])
ylim([-2,2])
axis('equal')
xlim([-2,2])
ylim([-2,2])

H_fh = matlabFunction(H);
figure, contour(xx1,xx2,H_fh(xx1,xx2),[0.05 0.6 0.7 0.8 0.9 0.99, 1.1])
hold on, plot([-1 1 1 -1],[1 1 -1 -1],'k.','MarkerSize',30)
axis equal
xlim([-2,2])
axis equal
xlim([-2,2])
ylim([-2,2])
