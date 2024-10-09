%% 
%  File: Lorenz_attractor.m
%  Directory: 7_ftools/workspace/07_DOA_positive_systems
%  Author: Peter Polcz (ppolcz@gmail.com) 
% 
%  Created on 2021. March 02. (2020b)

P_init(12.02)

%%
% Automatically generated stuff

G_reset(01)
%       │└─ verbosity (0:false, 1:true, 2:do not change)
%       └── scope depth (0:reset, 1:reset if not set)

He = he;
I = @(varargin) eye(varargin{:});
O = @(varargin) zeros(varargin{:});

%% model

x_lim_theoretical = [
    -1 1
    -1 1
    -1 1
    ];
[x_lfr,x_lfr_cell] = pcz_generateLFRStateVector('x',x_lim_theoretical);
pcz_generateSymStateVector(3,'x');
syms t real

rho = 28;
sigma = 10;
beta = 8/3;

Scale_x = 50;
Offset_x = [0,0,25]';

f = subs([
    sigma*(x2 - x1)
    x1*(rho - x3) - x2
    x1*x2 - beta*x3
    ],x,x*Scale_x+Offset_x)/Scale_x;

f_fh_cell = cellfun(@(f) {matlabFunction(f,'vars',x)}, num2cell(f));
f_lfr_cell = cellfun(@(f) {f(x_lfr_cell{:})},f_fh_cell);

% %%

f_ode = matlabFunction(f,'vars',{t, x});

%{

x0 = [1 1 1]';

fig = figure(2);
delete(fig.Children);
ax1 = axes('Parent',fig);
hold on

[tt,xx] = ode45(f_ode,[0 100],x0);
Pl1 = plot3(xx(:,1),xx(:,2),xx(:,3));
Pl2 = plot3(xx(1,1),xx(1,2),xx(1,3),'.', 'Color', Pl1.Color);

grid on
box on
axis equal

xlabel x
ylabel y
zlabel z

%}

%%

x1 = x_lfr_cell{1};
x2 = x_lfr_cell{2};
x3 = x_lfr_cell{3};
Pi_mon = [
    1
    x1
    x2
    x3
    ];
Pi_den = 1 + x1^2 + x2^2 + x3^2;
Pi_num = [
    1
    x1
    x2
    x3
    x1^2
    x2^2
    x3^2
    x1*x2
    x2*x3
    x3*x1
    % x1^3
    % x1^2*x2
    % x1^2*x3
    ];

Pi_lfr = [
    Pi_mon 
    Pi_num/Pi_den
    ];

% Pi_fh = matlabFunction(Pi,'vars',x,'Optimize',false);
% Pi_lfr = Pi_fh(x_lfr_cell{:});
Pi = plfr(Pi_lfr);
m = size(Pi_lfr,1);

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

x_lim = [
    -0.8 0.8
    -1 1
    -1 1
    ];

c1 = 0.1;
sdpvar c2;
% c2 = 0.9;
c3 = 1;

alpha0 = 3/8;

[X_v,~,X_fci,proj] = P_ndnorms_of_X(x_lim);

X_NrV = size(X_v,1);
[X_NrF,X_NrFc] = size(X_fci);

P = sdpvar(m);
L = sdpvar(m,s,'full');
Ld = sdpvar(md,sd,'full');

CONS = [ 0.00001 <= c2 , c2 <= c3*0.0199999];
OBJ = c2;

for i = 1:X_NrV
    xi = X_v(i,:)';
    
    CONS = [ CONS,
        P + He( L*N(xi) ) - 1e-7*eye(m) >= 0
        alpha0*He( Ed'*P*Ad + Ld*Nd(xi) ) + Ed'*P*Ed - blkdiag(c2,zeros(md-1)) + 1e-7*eye(md) <= 0
        ];
end

% Central condition
% xeq = [0.131,0.389,12.82]'/Scale_x;
% PiEq = Pi_plfr(xeq);
% CONS = [ CONS , PiEq' * P * PiEq <= c1 ];

% Boundary conditions
% tau = sdpvar(1,X_NrF);
Qf = @(alpha) P + blkdiag(-alpha,O(m-1));
for i = 1:X_NrF
    
    Verbosity = G_VERBOSE(0);
    [S_Fi,Pi_Fi,iS_Fi,Ker_Fi] = P_mingen_for_LFR(Pi,'proj',proj{i});
    N_Fi = P_affine_annihilator_for_LFR(Pi_Fi,x_lfr,'proj',proj{i});
    G_VERBOSE(Verbosity);
    
    L_b1 = sdpvar(N_Fi.nu,N_Fi.ny);
    % L_b2 = sdpvar(N_Fi.nu,N_Fi.ny);
    
    for j = 1:X_NrFc
        x_num = X_v(X_fci(i,j),:)';
        CONS = [ CONS
            S_Fi'*Qf(c3)*S_Fi + He( L_b1 * N_Fi(x_num) ) >= 0
            % S_Fi'*Qf(tau(i))*S_Fi + He( L_b2 * N_Fi(x_num) ) <= 0
            ];
    end
end

sol = optimize(CONS,OBJ,sdpsettings('solver','mosek'))
pcz_feasible(sol,CONS)

%%

c2 = double(c2);
L = double(L);
P = double(P);
Ld = double(Ld);

d = @(fun) jacobian(fun,x)*f;

phi0 = sym(plfr(Pi_mon));
dphi0 = d(phi0);

phi = sym(plfr(Pi_num));
dphi = d(phi);

q = sym(plfr(Pi_den));
dq = jacobian(q) * f;

m0 = numel(phi0);
m1 = numel(phi);
assert(size(P,1) == m0+m1);
Q = P(1:m0,1:m0);
R = P(1:m0,m0+1:m0+m1);
S = P(m0+1:m0+m1,m0+1:m0+m1);

p0 = phi0.' * Q * phi0;
p1 = phi0.' * R * phi;
p2 = phi.' * S * phi;

dp0 = d(p0);
dp1 = d(p1);
dp2 = d(p2);

p0_fh = matlabFunction(p0,'vars',x);
dp0_fh = matlabFunction(dp0,'vars',x);
p1_fh = matlabFunction(p1,'vars',x);
dp1_fh = matlabFunction(dp1,'vars',x);
p2_fh = matlabFunction(p2,'vars',x);
dp2_fh = matlabFunction(dp2,'vars',x);
q_fh = matlabFunction(q,'vars',x);
dq_fh = matlabFunction(dq,'vars',x);

resolution = 201;
xc = cellfun(@num2cell, num2cell(x_lim,2), 'UniformOutput', false);
lspace_cell = cellfun(@(lims) {linspace(lims{:},resolution)}, xc);

xx = cell(1,3);
[xx{:}] = ndgrid(lspace_cell{:});

pp0 = p0_fh(xx{:});
pp1 = p1_fh(xx{:});
pp2 = p2_fh(xx{:});
qq = q_fh(xx{:});

dpp0 = dp0_fh(xx{:});
dpp1 = dp1_fh(xx{:});
dpp2 = dp2_fh(xx{:});
dqq = dq_fh(xx{:});

VV = pp0 + 2*pp1 ./ qq + pp2 ./ (qq.^2);
dVV = dpp0 + 2*(dpp1.*qq - dqq.*pp1) ./ (qq.^2)  + dpp2 ./ (qq.^2) - 2*dqq.*pp2 ./ (qq.^3);

%%

Version = 1;

VV_half = VV;

switch Version
    case 1
        VV_half(xx{2} > 0 & xx{1} > 0) = NaN;
    case 2
        VV_half(xx{2} > 0 & xx{1} < 0) = NaN;
    case 3
        VV_half(xx{1} > xx{2}) = NaN;
end

Col = @(M) M(:);
V_level = min([
    Col(VV(:,:,[1 end]))
    Col(VV(:,[1,end],:))
    Col(VV([1,end],:,:))
    ]);

lambda = 0.00001;
c2_redef = (1-lambda)*max(VV(dVV(:) > 0)) + lambda*c2;

% %%

fig = figure(2);
delete(fig.Children);
ax1 = axes('Parent',fig);
hold on

% tmp = sqrt(beta*(rho-1));
% Color1 = pcz_get_plot_colors([],1);
% for x0 = ([0 0 0 ; tmp tmp rho-1 ; -tmp -tmp rho-1]'-Offset_x)/Scale_x
%     [~,x_ode] = ode89(f_ode,[0 100],x0);
%     Pl1 = plot3(x_ode(:,1),x_ode(:,2),x_ode(:,3),'Linewidth',2,'Color',Color1);
% end

Color1 = pcz_get_plot_colors([],1);
for x0 = [-0.1836 -0.0439 0.1971]'
    [~,x_ode] = ode89(f_ode,[0 100],x0);
    Pl1 = plot3(x_ode(:,1),x_ode(:,2),x_ode(:,3),'Color',Color1);
end

grid on
box on
axis equal

xlabel x
ylabel y
zlabel z

ph = patch(isosurface(xx{:},VV,c2_redef));
ph_Vertices = ph.Vertices;
ph_Faces = ph.Faces;
set(ph,'FaceColor','red','EdgeColor','none', 'FaceAlpha', 0.5);

ph = patch(isosurface(xx{:},VV_half,V_level));
ph_Vertices = ph.Vertices;
ph_Faces = ph.Faces;
set(ph,'FaceColor',[1 0.8431 0],'EdgeColor','none', 'FaceAlpha', 0.5);

ph = patch(isosurface(xx{:},dVV,0));
ph_Vertices = ph.Vertices;
ph_Faces = ph.Faces;
set(ph,'FaceColor','green','EdgeColor','none', 'FaceAlpha', 0.5);

% ph = patch(isosurface(xx{:},dVV_VV_1,0));
% ph_Vertices = ph.Vertices;
% ph_Faces = ph.Faces;
% set(ph,'FaceColor','blue','EdgeColor','none', 'FaceAlpha', 0.5);


axis equal
light,
switch Version
    case 1
        view([152.8652 25.4689])
    case 2
        view([-125.7518 9.8295])
    case 3
        view([-125.7518 9.8295])
end
xlabel $x_1$ interpreter latex
ylabel $x_2$ interpreter latex
zlabel $x_3$ interpreter latex
grid on
axlims = x_lim';
axis(axlims(:)),
box on
% set(gca,'LineWidth',1)

set(gca, Persist.font_axis12{:})
Persist.latexify_axis(gca,16)
Labels = Persist.latexified_labels(gca,20,'$x_1$','$x_2$','$x_3$');
Labels{1}.VerticalAlignment = 'bottom';
Labels{2}.VerticalAlignment = 'top';

good = VV < c2_redef;

VOLUME = prod(x_lim(:,2) - x_lim(:,1)) * ( sum(good(:))/numel(VV) )

% persist.print({'-dpdf'}, [ 'main_figure_pia1_' polytope '_without_legend.pdf' ])

Dirname = '/home/ppolcz/_/3_docs/28_nonzero_attractor/actual/fig/Exported_from_Matlab';
Figname = [ Dirname filesep datestr(date,29) '-Lorenz_attractor_LFR_min_c3_View1' ];

% pause(1)
% print([Figname '.png'],'-dpng','-r500')

return

%%

format long
for i = 1:X_NrV
    xi = X_v(i,:)';
    CHECK_LMI = alpha0*He( Ed'*P*Ad + Ld*Nd(xi) ) + Ed'*P*Ed - blkdiag(c2,zeros(md-1));

    eig(CHECK_LMI)
end
format short

xx_vec = [xx{1}(:) xx{2}(:) xx{3}(:)]; 
xx_invalid = xx_vec(dVV_VV_1(:) > 0,:);

xx_inv_cell = num2cell(num2cell(xx_invalid),2);

Pd = alpha0*He( Ed'*P*Ad ) + Ed'*P*Ed - blkdiag(c2,zeros(md-1));

Pid_fh(28.5,50,50)' * Pd_aug_fh(28.5,50,50) * Pid_fh(28.5,50,50)

check(CONS)

