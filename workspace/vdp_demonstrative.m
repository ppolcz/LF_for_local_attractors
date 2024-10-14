%%
%  File: ipend_rational_4D.m
%  Directory: 8_published/25_PhD_dissertation/1_DOA_computation
%  Author: Peter Polcz (ppolcz@gmail.com)
%
%  Created on 2020. May 19. (2019b)
%

P_init(12.01)
G_reset

He = he;
I = @(varargin) eye(varargin{:});

RUN_ID = pcz_runID(mfilename);
RUN_ID_str = num2str(RUN_ID);
setenv('RUN_ID', RUN_ID_str)
logger = Logger(['results/' mfilename '-output.txt']);
TMP_vcUXzzrUtfOumvfgWXDd = pcz_dispFunctionName;
pcz_dispFunction2('Run ID = %s', getenv('RUN_ID'));


rng(1);


%% Model description

% Generate symbolic state variables
P_generate_symvars(2,0,0,0);

x_lim_comp = [
    -1 1
    -1 1
    ]*10;
[xl,xl_cell] = pcz_generateLFRStateVector('x',x_lim_comp);

A_fh = @(x1,x2) [
    0          -1
    x1*x2+1    -1
    ];
    

A_sym = [
    0          -1
    prod(x)+1  -1
    ];
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
pcz_lfrzero_report(minlfr(S*pi - pi_init), 'Minimal generator for pi');

% Compute pi1 corresponding to the minimal generator pi
pi1 = iS(nx+1:end,nx+1:end) * pi1_init;
pcz_lfrzero_report(pi - [xl;pi1], 'Pi1 is');

m = pi.ny;

%% Structure of the Lyapunov function

% [DEFAULT]
phi = pi;
phi1 = pi1;
K = eye(m);
m_phi = phi.ny;

%{
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
pcz_lfrzero_report(minlfr(pi_new - pi), 'Find K: phi = K*pi, pi remained the same.')
pcz_lfrzero_report(minlfr(K*pi - phi), 'phi = K*pi')

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
pcz_lfrzero_report(pid_init - Sd*pid, 1e-6, 'Minimal generator for Pid')

% Annihilators
N = P_affine_annihilator_for_LFR(phi,xl);
Nd = P_affine_annihilator_for_LFR(pid,xl,'tol',1e-8);

pcz_lfrzero_report(N * phi,1e-10, 'N(x) phi(x) = 0')
pcz_lfrzero_report(Nd * pid, 1e-8, 'Nd(x) pid(x) = 0')

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

pcz_lfrzero_report(phi - Ed*pid, 1e-8, 'pi = Ed pid')
pcz_lfrzero_report(dphi - Ad*pid, 1e-6, 'dpi = Ad pid')

%%


X_v = [
   -0.2521    1.0387
   -1.4218    0.0202
   -0.5849   -1.2706
    1.0487   -0.8471
    1.3008    0.8370
    ];
[X_v,~,X_fci,proj] = P_2dnorms_of_X_v(X_v);

% Bounds
% x_lim = [
%     -1 1
%     -1 1
%     ] / 4;
% [X_v,~,X_fci,proj,~] = P_ndnorms_of_X(x_lim);
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

%     G_VERBOSE(0)
    pcz_dispFunction2('\nBoundary conditions %d/%d:',i,X_NrF);
    [S_Fi,Pi_Fi,iS_Fi,Ker_Fi] = P_mingen_for_LFR(phib,'proj',proj{i});
    N_Fi = P_affine_annihilator_for_LFR(Pi_Fi,xl,'proj',proj{i});
%     G_VERBOSE(1);

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

xorig_res = [51 51];
xorig_lim = [
    -1 1
    -1 1
    ]*2;

% Generate grid in X
lspace = cellfun(@(args) {linspace(args{:})}, num2cell(num2cell([xorig_lim xorig_res(:)]),2));
[x1_num,x2_num] = meshgrid(lspace{:});

% Vectorize grid in X
W_num = W_fh(x1_num,x2_num);

%%

fig = figure(1);
delete(fig.Children);
ax = axes('Parent',fig);

alpha = 0.9;

Min_Level = -alpha;
Max_Level = 3*alpha;

Sf = surf(x1_num,x2_num,W_num);
Sf.CData(W_num > Max_Level) = Max_Level;
Sf.CData = log(Sf.CData + 1);
colormap turbo

hold on
X_circ = [
    X_v
    X_v(1,:)
    ];

plot3(X_circ(:,1),X_circ(:,2),X_circ*[0;0] + Min_Level,'o-','LineWidth',2)

[Pc,fig] = pcontour(x1_num,x2_num,W_num,[1 1]*alpha);
delete(fig);

Fl = fill3(Pc.contour(1,:),Pc.contour(2,:),Pc.contour(2,:)*0 + Min_Level,...
    min(pcz_get_plot_colors([],2)*2,1));
Fl.FaceAlpha = 0.1;
Fl.EdgeAlpha = 0;
plot3(Pc.contour(1,:),Pc.contour(2,:),Pc.contour(2,:)*0 + Min_Level,...
    'LineWidth',2,'Color',pcz_get_plot_colors([],2));

xlim([-1.5 1.5])
ylim([-1.5 1.5])
zlim([Min_Level-eps Max_Level])

x0 = [
    0.8
    -0.4
    ];
f_ode = matlabFunction(f_sym,'vars',{t,x});
[tt,xx] = ode45(f_ode,[0,20],x0);

plot3(xx(:,1),xx(:,2),xx*[0;0] + Min_Level,'k')
plot3(xx(1,1),xx(1,2),Min_Level,'k.')
plot3(0,0,Min_Level,'k.','MarkerSize',15)

text(0,0,Min_Level,'origin','FontSize',12,'Interpreter','latex','VerticalAlignment','bottom','HorizontalAlignment','center')
text(xx(1,1),xx(1,2),Min_Level,'$x(0)$','FontSize',13,'Interpreter','latex','VerticalAlignment','top','HorizontalAlignment','right');
text(X_circ(4,1),X_circ(4,2),Min_Level,'$\mathcal X$','FontSize',15,'Interpreter','latex','VerticalAlignment','top','HorizontalAlignment','left')

box on

Logger.latexify_axis(gca,14)
Logger.latexified_labels(gca,14,'$x_1$','$x_2$','$V(x_1,x_2)$')

exportgraphics(gca,'/home/ppolcz/_/3_docs/30_ECC2021/fig/vdp_demonstrative.pdf','ContentType','vector')

return

%% Visualize 
% 
% [1] ASCII plot

% mat_fname ='/home/ppolcz/_/1_PhD_projects/25_PhD_dissertation/1_DOA_computation/results/ipend_rational_4D-output-2020-05-20_16:51_id2350_VEGLEGES/ipend4D_vel10_theta0.785_omega3_phi6.mat'; %#ok<UNRCH>
% mat_fname ='/home/ppolcz/_/1_PhD_projects/25_PhD_dissertation/1_DOA_computation/results/ipend_rational_4D-output-2020-05-26_13:57_id2378/ipend4D_vel15_theta0.785_omega3_phi15.mat'
% mat_fname ='/home/ppolcz/_/1_PhD_projects/25_PhD_dissertation/1_DOA_computation/results/ipend_rational_4D-output-2020-05-26_14:48_id2379/ipend4D_vel15_theta0.85_omega3_phi15.mat'
mat_fname ='/home/ppolcz/_/1_PhD_projects/25_PhD_dissertation/1_DOA_computation/results/ipend_rational_4D-output-2020-05-26_15:32_id2380/ipend4D_vel18_theta0.785_omega3_phi15.mat'
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

%% For LaTeX documentation

P_generate_symvars_v10(4,0,0,0)
syms v theta omega
x_expr = xorig2xz_fh(v,theta,omega).';
phi_inorig_sym = simplify(subs(phi_sym,x,x_expr))

pcz_latex_v7(phi_sym)

%%


% TMP_iwclITrbVyVrJaArrXNr = pcz_dispFunctionName('Evaluate the storage function over the grid');

% % Generate grid in X
% res = [ 31 ; 31 ; 31 ; 31 ];
% lspace = cellfun(@(args) {linspace(args{:})}, num2cell(num2cell([ x_lim res ]),2));
% x_grid = cell(1,nx);
% [x_grid{:}] = ndgrid(lspace{:});
% 
% V_num = V_plfr.val(x_grid{:});
% 
% pcz_dispFunctionEnd(TMP_iwclITrbVyVrJaArrXNr);
% 
% % -------------------------------------------------------------------------
% % -------------------------------------------------------------------------
% 
% TMP_GvDXGhRLfipwBoRPoGfI = pcz_dispFunctionName('Find maximal level set inside X');
% 
% % Find maximal level set inside X
% alpha = sdpvar;
% 
% % Hard coded (Now, I am a too lazy, to find a more adequate solution)
% CONS = [
%     V_num([1 end],:,:,:) >= alpha+eps
%     V_num(:,[1 end],:,:) >= alpha+eps
%     V_num(:,:,[1 end],:) >= alpha+eps
%     V_num(:,:,:,[1 end]) >= alpha+eps
%     ];
% 
% sol_alpha = optimize(CONS,-alpha,sdpsettings('verbose',0));
% alpha = double(alpha);
% 
% V_num01 = V_num;
% V_num01(V_num01 <= alpha) = 0;
% V_num01(V_num01 ~= 0) = 1;
% V_ascii = reshape(num2cell(char(V_num01*(' ' - '#') + '#'),[1 2]),[res(3)*res(4) 1]);
% 
% V_n0count = cellfun(@(v) sum(v(:)), num2cell(1-V_num01,[1 2]));
% V_n0count = V_n0count(:);
% z1_n0count = reshape(x_grid{3}(1,1,:,:),[res(3)*res(4) 1]);
% z2_n0count = reshape(x_grid{4}(1,1,:,:),[res(3)*res(4) 1]);
% 
% for i = 1:numel(V_n0count)
%     if V_n0count(i) < 1, continue, end
%     
%     fprintf('z1 = sin(theta) = %g, z2 = 1 - cos(theta) = %g\n', z1_n0count(i), z2_n0count(i))
%     disp(V_ascii{i})
% end
% 
% 
% pcz_dispFunction2('alpha = %g', alpha)
% 
% 
% pcz_dispFunctionEnd(TMP_GvDXGhRLfipwBoRPoGfI);
% 
% % -------------------------------------------------------------------------
% % -------------------------------------------------------------------------
