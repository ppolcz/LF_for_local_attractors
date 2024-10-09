%%
%  File: ipend_rational_4D.m
%  Directory: 8_published/25_PhD_dissertation/1_DOA_computation
%  Author: Peter Polcz (ppolcz@gmail.com)
%
%  Created on 2020. May 19. (2019b)
%

G_reset

He = he;
I = @(varargin) eye(varargin{:});

RUN_ID = pcz_runID(mfilename);
RUN_ID_str = num2str(RUN_ID);
setenv('RUN_ID', RUN_ID_str)
logger = Logger(['results/' mfilename '-output.txt']);
TMP_vcUXzzrUtfOumvfgWXDd = pcz_dispFunctionName;
pcz_dispFunction2('Run ID = %s', getenv('RUN_ID'));

% A few operations require random sampling. This line makes the results
% reproducible.
rng(1);


%% Model description

m = 1;       % Mass of the rod [kg]
M = 0.5;     % Mass of the car [kg]
L = 1;       % Length of the rod [m]
g = 10;      % Gravitational acceleration [m/s^2]
b = 0;       % Friction coefficient [kg/s]

% Moment of inertia of a rod of length 2L and mass m, rotating about one
% end. This expression assumes that the rod is an infinitely thin (but
% rigid) wire.
moment_of_In = 4*m*L^2/3;

% Generate symbolic state variables
P_generate_symvars(4,0,0,0);

% Named state variables
% y: not included
% v = x1;
% theta: not included
% omega = x2;
% z1 = x3; <--- sin(theta)
% z2 = x4; <--- 1 - cos(theta)

v_max = 18;
theta_max = builtin('pi')/4;
omega_max = 3;

x_lim = [
    -v_max          , v_max
    -omega_max      , omega_max
    -sin(theta_max) , sin(theta_max)
    -0.1            , 1-cos(theta_max) % this is technical (cannot be zero)
    ];

% Generate lfr state variables
[xl,xl_cell] = pcz_generateLFRStateVector('x',x_lim);

% Define auxiliary varibles
sigma1 = moment_of_In + m*L^2;
sigma3 = @(z1) sigma1*(m+M) - m^2*L^2*(1-z1^2);

K = [-1 -16.45 -38.56 0];

B_fh = @(v,omega,z1,z2) [
    sigma1
   -m*L*(1-z2)
    0
    0
    ] / sigma3(z1);

A_fh = @(v,omega,z1,z2) [
   -b*sigma1     ,  m*L*sigma1*omega*z1     , -m^2*L^2*g*(1-z2) , 0
    b*m*L*(1-z2) , -m^2*L^2*omega*z1*(1-z2) ,  (m+M)*m*g*L      , 0
    0            ,  sigma3(z1)*(1-z2)       ,  0                , 0
    0            ,  0                       ,  sigma3(z1)*omega , 0
    ] / sigma3(z1) - B_fh(v,omega,z1,z2) * K;

A_sym = simplify(A_fh(x_cell{:}));
f_sym = simplify(A_sym * x);

%%

fi_fh_cell = cellfun(@(fi) {matlabFunction(fi,'vars',x,'Optimize',false)},num2cell(f_sym));
fi_lfr_cell = cellfun(@(fi) {fi(xl_cell{:})}, fi_fh_cell);

%% Model transformation

A_plfr = plfr(A_fh(xl_cell{:}));
A_plfr = A_plfr.set_vars(x);

% Initial generator pi1 and pi
pi_init = A_plfr.generatePI * xl;
pi1_init = A_plfr.generatePI1 * xl;

% Minimal generator for pi
[S,pi,iS,Ker] = P_mingen_for_LFR(pi_init);
warning pcz_gsszero_report
pcz_gsszero_report(minlfr(S*pi - pi_init), 'Minimal generator for pi');

% Compute pi1 corresponding to the minimal generator pi
pi1 = iS(nx+1:end,nx+1:end) * pi1_init;
pcz_gsszero_report(pi - [xl;pi1], 'Pi1 is');

m = pi.ny;

%% Structure of the Lyapunov function

% [DEFAULT]
phi = pi;
phi1 = pi1;
K = eye(m);
m_phi = phi.ny;

% %{
phi1_fh = @(x1,x2,x3,x4) [
      x1*x3 % deg: 2
      x1*x4 % deg: 2
      x2*x3 % deg: 2
      x2*x4 % deg: 2
      x3*x4 % deg: 2
       x3^2 % deg: 2
       ...
    x1*x3^2 % deg: 3
   x2*x3*x4 % deg: 3
    x2*x3^2 % deg: 3
    x2^2*x3 % deg: 3
    x3^2*x4 % deg: 3
       x3^3 % deg: 3
       ...
    x2*x3^3 % deg: 4
 x2*x3^2*x4 % deg: 4
 x2^2*x3*x4 % deg: 4
    ] / (2*x3^2 + 5);

phi1 = plfr(phi1_fh(xl_cell{:}));

phi = [
    plfr(xl)
    phi1
    ];
m_phi = phi.ny;

% Find K:
[K_all,pi_new,iK,Ker] = P_mingen_for_LFR([pi;phi]);
K = K_all(m+1:m+m_phi,:);
pcz_gsszero_report(minlfr(pi_new - pi), 'Find K: phi = K*pi, pi remained the same.')
pcz_gsszero_report(minlfr(K*pi - phi), 'phi = K*pi')

%}

%%

% Compute the time derivative of phi1 and phi
dphi1 = diff(phi1,xl_cell,fi_lfr_cell);
dphi = [ A_plfr * xl ; dphi1 ];

% Initial generator Pid
pid_init = [
    pi
    dphi1
    ];

% Inflate Pid using the LFR approach (Method II of [Polcz etal, EJC 2018])
PI_of_pid_init = pid_init.generatePI;
F12_of_pid_init = [ pid_init.A pid_init.B ];
[Sd_init,pid,iSd,Ker_d] = P_mingen_for_LFR(PI_of_pid_init);
Sd = F12_of_pid_init * Sd_init;
pcz_gsszero_report(pid_init - Sd*pid, 1e-6, 'Minimal generator for Pid')

% Annihilators
N = P_affine_annihilator_for_LFR(phi,xl);
Nd = P_affine_annihilator_for_LFR(pid,xl,'tol',1e-8);

pcz_gsszero_report(N * phi,1e-10, 'N(x) phi(x) = 0')
pcz_gsszero_report(Nd * pid, 1e-8, 'Nd(x) pid(x) = 0')

% Sizes
m1_phi = phi1.ny;
m1 = pi1.ny;
s = N.ny;
md = pid.ny;
sd = Nd.ny;

Ed = [ K zeros(m_phi,m1_phi) ] * Sd;
Ad = [
    [ A_plfr.A A_plfr.B ] * S , zeros(nx,m1_phi)
    zeros(m1_phi,m)           , eye(m1_phi)
    ] * Sd;

pcz_gsszero_report(phi - Ed*pid, 1e-8, 'pi = Ed pid')
pcz_gsszero_report(dphi - Ad*pid, 1e-6, 'dpi = Ad pid')

%%

% Bounds
[X_v,X_ak,X_fci,proj,~] = P_ndnorms_of_X(x_lim);
X_NrV = size(X_v,1);
[X_NrF,X_NrFc] = size(X_fci);

P = sdpvar(m_phi);
L = sdpvar(m_phi,s,'full');
Ld = sdpvar(md,sd,'full');

CONS = [];

CONS_desc_str = {
    '        V > 0       '
    '       dV > 0       '
    ' boundary condition '
    ' b.c. with slack v. '
    };
CONS_desc = {
    '++++++++++++++++++++'
    '   Condition type   '
    '++++++++++++++++++++'
    };

for i = 1:X_NrV
    xi = X_v(i,:)';

    CONS = [ CONS,
        P + He( L*N(xi) ) - 1e-5*I(m_phi) >= 0
        He( Ed'*P*Ad + Ld*Nd(xi) ) + 0*I(md) <= 0
        ];
    CONS_desc = [ CONS_desc ; CONS_desc_str(1:2)];
end

% Boundary conditions
tau = sdpvar(1,X_NrF);
Qf = @(alpha) blkdiag(-alpha,P);
phib = [1;phi];
for i = 1:X_NrF

    % G_VERBOSE(0)
    pcz_dispFunction2('\nBoundary conditions %d/%d:',i,X_NrF);
    [S_Fi,Pi_Fi,iS_Fi,Ker_Fi] = P_mingen_for_LFR(phib,'proj',proj{i});
    N_Fi = P_affine_annihilator_for_LFR(Pi_Fi,xl,'proj',proj{i});
    % G_VERBOSE(1);

    L_b1 = sdpvar(N_Fi.nu,N_Fi.ny);
    L_b2 = sdpvar(N_Fi.nu,N_Fi.ny);

%{
    % Ellenorzeskeppen
    Ker_Fi = round(rref(Ker_Fi),10)
    N_Fi = P_affine_annihilator_for_LFR(Pi_Fi,xl,'proj',proj{i},'sym',1);
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
        
        CONS_desc = [ CONS_desc ; CONS_desc_str(3:4)];
            
    end
end
CONS_desc = [ CONS_desc ; CONS_desc(1) ];

OBJ = sum(tau);

sol = optimize(CONS,OBJ)
pcz_feasible(sol,CONS)

%% Save results

TMP_SBWFrkjLRJdKsBtETTdR = pcz_dispFunctionName('Save results');


disp_CONS = strsplit(evalc('display(CONS)'),newline)';
check_CONS = strsplit(evalc('check(CONS)'),newline)';

cnum = numel(disp_CONS)-1;

str_CONS = num2cell([ disp_CONS(1:cnum) CONS_desc(1:cnum) check_CONS(2:cnum+1) ],2);
str_CONS = cellfun(@(c) { horzcat(c{1}(1:end),c{2},c{3}(30:end)) },str_CONS);
str_CONS = strjoin(str_CONS,newline);
pcz_dispFunction2(evalc('disp(str_CONS)'));

P_val = value(P);
V_plfr = phi'* P_val * phi;
V_plfr_M = V_plfr.M;
V_plfr_Delta = V_plfr.Delta;
V_plfr_blk = V_plfr.blk;
phi_sym = sym(phi);
pcz_dispFunction2(evalc('display(phi_sym)'))

W_sym = phi_sym.' * P_val * phi_sym;
W_fh = matlabFunction(W_sym,'vars',x);

dW_sym = jacobian(W_sym,x) * f_sym;
dW_fh = matlabFunction(dW_sym,'vars',x);

xorig_res = [71 71 71];
xorig_lim = [
    -v_max     , v_max
    -theta_max , theta_max
    -omega_max , omega_max
    ];

% Generate grid in X
lspace = cellfun(@(args) {linspace(args{:})}, num2cell(num2cell([xorig_lim xorig_res(:)]),2));
xorig_grid = cell(1,size(xorig_lim,1));
[xorig_grid{:}] = ndgrid(lspace{:});

% Vectorize grid in X
grid_size = size(xorig_grid{1});
xorig_col = cellfun(@(a) {a(:)}, xorig_grid);

% Compute the values of p over the grid in X
xorig2xz_fh = @(v,theta,omega) [ v omega sin(theta) 1-cos(theta) ];
xz_col = num2cell(xorig2xz_fh(xorig_col{:}),1);

% Compute the values of the storage function
W_evaluation_Timer = pcz_dispFunctionName('Evaluate the storage function over the grid');
W_col = W_fh(xz_col{:});
V_num = reshape(W_col,grid_size);
pcz_dispFunctionEnd(W_evaluation_Timer);

% Compute the values of the derivative function
dW_evaluation_Timer = pcz_dispFunctionName('Evaluate the derivative function over the grid');
dW_col = dW_fh(xz_col{:});
dV_num = reshape(dW_col,grid_size);
pcz_dispFunctionEnd(dW_evaluation_Timer);

alpha = 1;
VOL_X = prod(xorig_lim(:,2) - xorig_lim(:,1));
VOL_Gamma = VOL_X * ( sum(V_num(:) < alpha) / numel(V_num) );
pcz_dispFunction2('Volume of the X: %g', VOL_X)
pcz_dispFunction2('Approximated volume of the invariant domain: %g', VOL_Gamma)

V_num01 = V_num;
V_num01(V_num01 <= alpha) = 0;
V_num01(V_num01 ~= 0) = 1;
V_ascii = reshape(num2cell(char(V_num01*(' ' - '#') + '#'),[1 2]),[xorig_res(3) 1]);

dV_num01 = dV_num;
dV_num01(dV_num01 <= 0) = 0;
dV_num01(dV_num01 > 0) = 1;
dV_ascii = reshape(num2cell(char(dV_num01*(' ' - '#') + '#'),[1 2]),[xorig_res(3) 1]);

V_n0count = cellfun(@(v) sum(v(:)), num2cell(1-V_num01,[1 2]));
V_n0count = V_n0count(:);
theta_n0count = reshape(xorig_grid{2}(1,:,1),[xorig_res(3) 1]);

for i = 1:numel(V_n0count)
    if V_n0count(i) < 1, continue, end
    
    pcz_dispFunction2('\ntheta = %g, velocity(%g↑%d) omega(%g→%g) :\n', theta_n0count(i),xorig_lim(1,:),xorig_lim(3,:))
    V_ascii{i}([1,end],:) = '─';
    V_ascii{i}(:,[1,end]) = '│';
    V_ascii{i}(1,1) = '└';
    V_ascii{i}(1,end) = '┘';
    V_ascii{i}(end,end) = '┐';
    V_ascii{i}(end,1) = '┌';
    V_ascii_cell = num2cell(V_ascii{i},2);
    V_ascii_cell{1} = [V_ascii_cell{1} ' ' num2str(xorig_lim(1,1))];
    V_ascii_cell{end} = [V_ascii_cell{end} ' ' num2str(xorig_lim(1,2))];
    V_ascii_cell(:) = V_ascii_cell(end:-1:1);
    V_ascii{i} = strjoin(V_ascii_cell,newline);
    pcz_dispFunction2(V_ascii{i})
end


% Save results (.mat)

w_ticks = -omega_max:1:omega_max;
v_ticks = -v_max:5:v_max;
theta_ticks = [-theta_max theta_max];
theta_ticklabels = {'$-\bar{\vartheta}$','$\bar{\vartheta}$'};

modelname = sprintf('ipend4D_vel%g_theta%g_omega%g_phi%d',v_max,round(theta_max,3),omega_max,size(phi1,1));
mat_fname = logger.mat_fname(modelname);
pcz_dispFunction2(mat_fname)
pcz_save(mat_fname, modelname,v_max,theta_max,omega_max,xorig_res,xorig_lim,xorig_grid,lspace,grid_size,...
    xorig2xz_fh,xz_col,alpha,V_num01,V_ascii,V_n0count,V_num,theta_n0count,VOL_X,VOL_Gamma,P_val,...
    w_ticks,v_ticks,theta_ticks,theta_ticklabels,P_val,V_plfr_M,V_plfr_Delta,V_plfr_blk,...
    phi_sym,logger);

pcz_dispFunctionEnd(TMP_SBWFrkjLRJdKsBtETTdR);

pcz_dispFunctionEnd(TMP_vcUXzzrUtfOumvfgWXDd);
logger.stoplog
return

%% Visualize 
% 
% [1] ASCII plot



mat_fname ='results/ipend_rational_4D-output-2020-09-30_18:37_id0005/ipend4D_vel18_theta0.785_omega3_phi15.mat';
load(mat_fname)

for i = 1:numel(V_n0count)
    if V_n0count(i) < 1, continue, end    
    pcz_dispFunction2('\ntheta = %g, velocity(%g↑%d) omega(%g→%g) :\n', theta_n0count(i),xorig_lim(1,:),xorig_lim(3,:))
    pcz_dispFunction2(V_ascii{i})
end


% [2] Nice viusal plot (just load)

TMP_oyXrvtrEzBjWNTyEyMog = pcz_dispFunctionName('Generate iso-surface of the Lyapunov/storage function');

fig = gcf;
delete(fig.get('Children'));

fontsize = 14;
v_ticks = -v_max:6:v_max;

...........................................................................

ax1 = subplot(1,5,1:4);
ph = patch(isosurface(xorig_grid{:},V_num,alpha));
ph_Vertices = ph.Vertices;
ph_Faces = ph.Faces;
set(ph,'FaceColor','red','EdgeColor','none', 'FaceAlpha', 0.8);

axis equal
light('Position',[30 15 -5]),
view([ -135.380 , 13.129 ])

grid on
axis(reshape(xorig_lim',[numel(xorig_lim) 1])),
box on

set(gca, Logger.font_axis12{:})
Logger.latexify_axis(gca,fontsize)
Labels_ax1 = Logger.latexified_labels(gca,fontsize,'$v$','','$\omega$');
Labels_ax1{1}.Position = [0 0 -3.7665];
xticks(v_ticks);
yticks(theta_ticks);
yticklabels(theta_ticklabels)
zticks(w_ticks);

rotate3d on

...........................................................................

ax1 = subplot(1,5,5);
ph = patch(isosurface(xorig_grid{:},V_num,alpha));
ph_Vertices = ph.Vertices;
ph_Faces = ph.Faces;
set(ph,'FaceColor','red','EdgeColor','none', 'FaceAlpha', 0.8);


axis equal
light('Position',-[-30 -15 5]),
view([-90 0])

grid on
axis(reshape(xorig_lim',[numel(xorig_lim) 1])),
box on
% set(gca,'LineWidth',1)

set(gca, Logger.font_axis12{:})
Logger.latexify_axis(gca,fontsize)
Labels_ax2 = Logger.latexified_labels(gca,fontsize,'$v$','','$\omega$');
xticks(v_ticks);
yticks(theta_ticks);
yticklabels(theta_ticklabels)
zticks(w_ticks);

rotate3d on

pcz_dispFunctionEnd(TMP_oyXrvtrEzBjWNTyEyMog);

%{

print(logger.fig_fname('main.png'),'-dpng','-r500')

%}

