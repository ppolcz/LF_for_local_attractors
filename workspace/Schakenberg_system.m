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
% verbosity (0:false, 1:true, 2:do not change)
% scope depth (0:reset, 1:reset if not set)

He = he;
I = @(varargin) eye(varargin{:});
O = @(varargin) zeros(varargin{:});

%% model

x_lim_theoretical = [
    -1 1
    -1 1
    ];
[x_lfr,x_lfr_cell] = pcz_generateLFRStateVector('x',x_lim_theoretical);
pcz_generateSymStateVector(2,'x');
syms t real

a = 1/8;
mu = 0.5;

f = [
    a - x1 + x1^2*x2
    mu - x1^2*x2
    ];
f_ode = matlabFunction(f,'Vars',{t,x});
[~,xx] = ode89(f_ode,[0 100],[1;1]);
xEnd = xx(end,:)';

%%

Pi = pcz_monomials(x,0:3);
Pi_fh = matlabFunction(Pi,'vars',x,'Optimize',false);
Pi_lfr = Pi_fh(x_lfr_cell{:});
Pi_plfr = plfr(Pi_lfr);
m = numel(Pi);

Pid = pcz_monomials(x,0:5);
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

%%

x_lim = [
    -4 6
    -4 6
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

CONS = [ c1*1.00001 <= c2 , c2 <= c3*0.99999];
OBJ = [];

for i = 1:X_NrV
    xi = X_v(i,:)';
    
    CONS = [ CONS,
        P + He( L*N(xi) ) - 1e-7*eye(m) >= 0
        alpha0*He( Ed'*P*Ad + Ld*Nd(xi) ) + Ed'*P*Ed - blkdiag(c2,zeros(md-1)) + 1e-7*eye(md) <= 0
        ];
end

% Central condition
PiEq = Pi_plfr(xEnd);
CONS = [ CONS , PiEq' * P * PiEq <= c1 ];

% Boundary conditions
% tau = sdpvar(1,X_NrF);
Qf = @(alpha) P + blkdiag(-alpha,O(m-1));
for i = 1:X_NrF
    
    Verbosity = G_VERBOSE(0);
    [S_Fi,Pi_Fi,iS_Fi,Ker_Fi] = P_mingen_for_LFR(Pi_plfr,'proj',proj{i});
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
Pi_sym = sym(Pi);
Pid_sym = sym(Pid);

V = Pi_sym'*P*Pi_sym;
dV = Pid' * He( Ed'*P*Ad ) * Pid;
Pd = alpha0*He( Ed'*P*Ad ) + Ed'*P*Ed - blkdiag(c2,zeros(md-1));
Pd_aug = alpha0*He( Ed'*P*Ad + Ld*sym(Nd) ) + Ed'*P*Ed - blkdiag(c2,zeros(md-1));
dV_V_1 = Pid' * Pd * Pid;

V_fh = matlabFunction(V);
dV_fh = matlabFunction(dV);
dV_V_1_fh = matlabFunction(dV_V_1);
Pid_fh = matlabFunction(Pid);
Pd_aug_fh = matlabFunction(Pd_aug);

resolution = 201;
xc = cellfun(@num2cell, num2cell(x_lim,2), 'UniformOutput', false);
lspace_cell = cellfun(@(lims) {linspace(lims{:},resolution)}, xc);

xx = cell(1,2);
[xx{:}] = ndgrid(lspace_cell{:});

VV = V_fh(xx{:});
dVV = dV_fh(xx{:});
dVV_VV_1 = dV_V_1_fh(xx{:});

VV_half = VV;
VV_half(xx{2} > 0 & xx{1} > 0) = NaN;

Col = @(M) M(:);
V_level = min([
    Col(VV(:,[1,end]))
    Col(VV([1,end],:))
    ]);

lambda = 0.01;
c2_redef = (1-lambda)*max(VV(dVV(:) > 0)) + lambda*c2;

%%

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
for x0 = [1;1]
    [~,x_ode] = ode89(f_ode,[0 100],x0);
    Pl1 = plot(x_ode(:,1),x_ode(:,2),'Color',Color1);
end

grid on
box on

xlabel x
ylabel y
zlabel z

TrV = @(V) log10(V + 1 - min(VV(:))) + 1;

surf(xx{:},TrV(VV))
contour(xx{:},TrV(VV),(1-logspace(-4,-1.5,20))*TrV(V_level))
view([65 20])
shading interp
axis tight
axis vis3d

