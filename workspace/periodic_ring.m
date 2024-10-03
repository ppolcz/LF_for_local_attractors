%% 
%  
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  Created on 2021. March 25. (2020b)
%  Minor review on 2021. June 17. (2021a)
%
% EXECUTIONS:
% - 2024.10.03. (October  3, Thursday), 11:14 (Matlab R2023b)

%%
% Automatically generated stuff

G_reset(01)
%        verbosity (0:false, 1:true, 2:do not change)
%        scope depth (0:reset, 1:reset if not set)

He = he;
I = @(varargin) eye(varargin{:});
O = @(varargin) zeros(varargin{:});

P_init(12)

%% model

alpha0 = 0.1;
c1 = 0.1;
c2 = 1;
c3 = 1;

% 2021.02.25. (februr 25, cstrtk), 14:29
Radius = 2;
x_lim = [
    -1 1
    -1 1
    ] * 3;
[x_lfr,x_lfr_cell] = pcz_generateLFRStateVector('x',x_lim);
pcz_generateSymStateVector(2,'x');
syms t real

mu = 1;

f = [
    x2 - mu*x1*(x1^2 + x2^2 - Radius^2)
    -x1 - mu*x2*(x1^2 + x2^2 - Radius^2)
    ];

%%

f_ode = matlabFunction(f,'vars',{t, x});

%{
x0 = [-0.2 0.39246321]';

fig = figure(2);
delete(fig.Children);
ax1 = axes('Parent',fig);
hold on

[tt,xx] = ode45(f_ode,[0 100],x0);
Pl1 = plot(xx(:,1),xx(:,2));
Pl2 = plot(xx(1,1),xx(1,2),'.', 'Color', Pl1.Color);

axis equal
%}
%%

Pi = pcz_monomials(x,0:2);
Pi_fh = matlabFunction(Pi,'vars',x,'Optimize',false);
Pi_lfr = Pi_fh(x_lfr_cell{:});
Pi_plfr = plfr(Pi_lfr);
m = numel(Pi);

Pid = pcz_monomials(x,0:4);
Pid_fh = matlabFunction(Pid,'vars',x,'Optimize',false);
Pid_lfr = Pid_fh(x_lfr_cell{:});
Pid_plfr = plfr(Pid_lfr);
md = numel(Pid);

N = P_affine_annihilator(Pi,x,'sym',1);
Nd = P_affine_annihilator(Pid,x,'sym',1);

% Compute Ad: dPi = Ad*Pid
dPi = jacobian(Pi,x)*f;
Fd = sym('F',[m md]);
[c,~] = pcoeffs(dPi - Fd*Pid,x);
[AA,bb] = equationsToMatrix(horzcat(c{:}),Fd(:));
Ad = double(reshape(AA\bb,size(Fd)));

pcz_symzero(dPi - Ad*Pid, 'dot pi(x) = Ad pi(x)');

% Pi = Ed*Pid
Ed = [ eye(m) zeros(m,md-m) ];

[s,m] = size(N);
[sd,md] = size(Nd);

% %{
%% Ez mar konkretan a cikkhez kell

pcz_latex_v7(sym(N))
pcz_latex_v7(sym(Nd))

%%
%}


%%

[X_v,~,X_fci,proj] = P_ndnorms_of_X(x_lim);



% X_v = [
%    -2.5    1.0
%     1.0    3.0
%     2.5    1.0
%     2.5   -1.0
%    -1.0   -3.0
%    -2.5   -1.0
%    ];
% [X_v,~,X_fci,proj] = P_ndnorms_of_X([],'Vertices',X_v);


X_NrV = size(X_v,1);
[X_NrF,X_NrFc] = size(X_fci);

P = sdpvar(m);
L = sdpvar(m,s,'full');
Ld = sdpvar(md,sd,'full');

alpha0 = 0.1;

lambda = 0.1;
sdpvar tau0
CONS = tau0 <= ((1-lambda)*c1+lambda*c2)/alpha0;

for i = 1:X_NrV
    xi = X_v(i,:)';
    
    CONS = [ CONS,
        P + He( L*N(xi) ) - 1e-5*I(m) >= 0
        He( Ed'*P*Ad + Ld*Nd(xi) ) + 0*I(md) <= 0
        % He( Ed'*P*Ad + Ld*Nd(xi) ) + alpha0*Ed'*P_store*Ed - blkdiag(tau0,zeros(md-1)) <= 0
        ];
end

% Central condition
xeq = [2 0]';
PiEq = Pi_plfr(xeq);
CONS = [ CONS , PiEq' * P * PiEq <= c1 ];


% Boundary conditions
tau = sdpvar(1,X_NrF);
Qf = @(alpha) P + blkdiag(-alpha,O(m-1));
Qfd = @(alpha) He( Ed'*P*Ad ) + blkdiag(-alpha,O(md-1));
for i = 1:X_NrF
    
    Verbosity = G_VERBOSE(0);
    [S_Fi,Pi_Fi,iS_Fi,Ker_Fi] = P_mingen_for_LFR(Pi_plfr,'proj',proj{i});
    N_Fi = P_affine_annihilator_for_LFR(Pi_Fi,x_lfr,'proj',proj{i},'sym',1);

    [S_Fid,Pi_Fid,iS_Fid,Ker_Fid] = P_mingen_for_LFR(Pid_plfr,'proj',proj{i});
    N_Fid = P_affine_annihilator_for_LFR(Pi_Fid,x_lfr,'proj',proj{i},'sym',1);
    G_VERBOSE(Verbosity);
    
    %{
    %% Ez mar konkretan a cikkhez kell

    pcz_num2str_latex(S_Fi);
    pcz_latex_v7(sym(Pi_Fi),'label',sprintf('\\\\pi_{f,%d} = ',i))
    pcz_latex_v7(sym(N_Fi),'label',sprintf('N_{f,%d} = ',i))

    pcz_num2str_latex(S_Fid);
    pcz_latex_v7(sym(Pi_Fid),'label',sprintf('\\\\pi_{f,%d} = ',i))
    pcz_latex_v7(sym(N_Fid),'label',sprintf('N_{f,%d} = ',i))

    %%
    %}

    L_b1 = sdpvar(N_Fi.nu,N_Fi.ny);
    L_b2 = sdpvar(N_Fi.nu,N_Fi.ny);

    L_b1d = sdpvar(N_Fid.nu,N_Fid.ny);
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
            S_Fi'*Qf(c2)*S_Fi + He( L_b1 * N_Fi(x_num) ) >= 0
            % S_Fi'*Qf(tau(i))*S_Fi + He( L_b2 * N_Fi(x_num) ) <= 0

            % S_Fid'*Qfd(-c3)*S_Fid + He( L_b1d * N_Fid(x_num) ) <= 0            
            ];
    end
end

% OBJ = sum(tau);
OBJ = [];
% OBJ = tau0;

sol = optimize(CONS,OBJ)
pcz_feasible(sol,CONS)

%% Visualization

L = double(L);
Ld = double(Ld);
P = double(P);
Pi_sym = sym(Pi);
Pid_sym = sym(Pid);

%{
%% Ez mar konkretan a cikkhez kell

pcz_num2str_latex(P);

%%
%}

[xx1,xx2] = meshgrid(linspace(x_lim(1,1),x_lim(1,2),41),linspace(x_lim(2,1),x_lim(2,2),41));

dV = Pid' * He( Ed'*P*Ad ) * Pid;

% 2021.02.26. (februr 26, pntek), 09:03
% p = Pi_sym'*P_store*Pi_sym;
% p_fh = matlabFunction(p);
% p_num = p_fh(xx1,xx2);
% p_level = value(tau0) / alpha0;

% 2021.02.26. (februr 26, pntek), 09:15
V = Pi_sym'*P*Pi_sym;
V_fh = matlabFunction(V);
V_num = V_fh(xx1,xx2);
V_level = min(min([
    V_num(:,[1 end])
    V_num([1,end],:)'
    ]));

% dV_fh = matlabFunction(jacobian(V)*f);
dV_fh = matlabFunction(dV);
dV_num = dV_fh(xx1,xx2);

resolution = { 471 471 }.';

lscell = cellfun(@(c) {linspace(c{:})}, num2cell([ num2cell(x_lim) resolution ],2) );
xxV = cell(1,3);
[xxV{1:2}] = ndgrid(lscell{:});
xxdV = xxV;
xxp = xxV;

% Retain invariant level set
figtmp = figure;
% xxp{3} = p_fh(xxp{1:2});
% Cont_p = pcontour(xxp{:},[1 1]*p_level);
xxV{3} = V_fh(xxV{1:2});
Cont_V = pcontour(xxV{:},[1 1]*V_level);
xxdV{3} = dV_fh(xxdV{1:2});
Cont_dV = pcontour(xxdV{:},[0 0]);
close(figtmp);

figtmp = figure;
xxV{3} = V_fh(xxV{1:2});
Cont_V_multiple = pcontour(xxV{:},[0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9]);
close(figtmp);

%% Fancy plot a Lyapunov fuggveny

fig = figure(6);
delete(fig.Children);
ax1 = subplot(121);

Level_2D = -5;
Threshold_V = 4;
Threshold_dV = -100;

V_clip = V_num;
V_clip(V_clip > Threshold_V*2) = nan;
V_clip(xx1 > 0 & xx2 < 0) = nan;

dV_clip = dV_num;
% dV_clip(dV_clip < Threshold_dV*2) = nan;

figure(fig);
hold off;
Surf = surf(xx1,xx2,V_clip);
Surf.CData = log10(V_clip);
Surf.FaceColor = 'interp';
hold on

% Limit cycle
[tt,xx] = ode45(f_ode,0:0.01:7,[Radius 0]');
Pl1 = plot3(xx(:,1),xx(:,2),xx(:,2)*0+Level_2D,'k','LineWidth',3);

V_along_limit_cycle = V_fh(xx(:,1),xx(:,2));

[tt,xx] = ode45(f_ode,[0 7],[0 0]');
% Pl2 = plot3(xx(:,1),xx(:,2),xx(:,2)*0+Level_2D,'.k','LineWidth',3,'MarkerSize',30);

for C = Cont_V
    c = C.contour;
    plot3(c(1,:),c(2,:),c(2,:)*0+Level_2D,'LineWidth',2,'Color',pcz_get_plot_colors([],1))
end

TSurf = fill3(c(1,:),c(2,:),c(2,:)*0+Level_2D,pcz_get_plot_colors([],1));
TSurf.FaceAlpha = 0.3;

axis(reshape(x_lim',[1 4]))
grid on

ax1.View = [ 73.4658   48.8985 ];
ax1.XTick = [x_lim(1,1):x_lim(1,2)];
ax1.YTick = [x_lim(2,1):x_lim(2,2)];
ax1.ZTick = 0:2:Threshold_V;
ax1.ZLim = [Level_2D Threshold_V];

axis vis3d
colormap jet


X_ind = convhull(X_v(:,1),X_v(:,2));
X_v1 = X_v(X_ind,:);
X1 = plot3(X_v1(:,1),X_v1(:,2),X_v1(:,1)*0+Level_2D,'o-','Color',pcz_get_plot_colors([],2),'LineWidth',1.5);

Logger.latexified_labels(ax1,22,'$x_1$','$x_2$');
PersistStandalone.latexify_axis(18)

Text3_V = text(-4.5,-2.5,Threshold_V,'$V(x)$',PersistStandalone.font_latex18{:});
Text3_Omega_cM = text(0,-1,Level_2D+1,'$\Omega_{c,M}$',PersistStandalone.font_latex18{:});
Text3_M = text(1.9,-1,Level_2D+1,'$M$',PersistStandalone.font_latex18{:});

ax2 = subplot(122);
Surf_dV = surf(xx1,xx2,dV_clip);
Surf_dV.CData = -log(abs(dV_clip)+1) .* sign(dV_clip);
Surf_dV.FaceColor = 'interp';

axis(reshape(x_lim',[1 4]))
ax2.XTick = [x_lim(1,1):x_lim(1,2)];
ax2.YTick = [x_lim(2,1):x_lim(2,2)];
ax2.ZLim = [Threshold_dV max(dV_clip(:))];

ax2.View = [ 73.4658   48.8985 ];

axis vis3d

Logger.latexified_labels(ax2,22,'$x_1$','$x_2$');
Logger.latexify_axis(18)


Text3_dV = text(-4.5,-2.5,10,'$\frac{\partial V}{\partial x}(x) f(x)$',PersistStandalone.font_latex18{:});

%{

Dirname = '/home/ppolcz/T/Dropbox/Peti/Munka/01_PPKE_2020/Dokumentaciok/28_nonzero_attractor_actual/fig/Exported_from_Matlab';
Figname = [ Dirname filesep datestr(date,29) '-periodic_ring_4th_ord' ];

pause(1)
print([Figname '.png'],'-dpng','-r500')

pngname = persist.file_simple('fig','vdp_10_infeasible_4ed_foku.png');
print(pngname,'-dpng','-r500')

%}

