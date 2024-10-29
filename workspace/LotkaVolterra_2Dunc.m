%%
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  Created on 2024. October 24. (2023a)
%

P_init(12.02)

%%
% Automatically generated stuff

G_reset(01)
%       │└─ verbosity (0:false, 1:true, 2:do not change)
%       └── scope depth (0:reset, 1:reset if not set)

He = he;
I = @(varargin) eye(varargin{:});
O = @(varargin) zeros(varargin{:});

%%

[xg,xg_cell] = pcz_generateGSStateVector('x',[0.1 4 ; 0.1 4],[-5,5 ; -5,5]);
[pg,pg_cell] = pcz_generateGSStateVector('p',[1.8 2.2 ; 0.8 1.0],[-5,5 ; -5,5]);

a1 = pg_cell{1};
a4 = pg_cell{2};

A_lv = [
    -a1 -3
    1.4 a4
    ];

b_lv = [
    5
    -2.4
    ];

xeq = -A_lv \ b_lv;

xhat = xg + xeq;

f = plfr([xhat(1) 0 ; 0 xhat(2)]*A_lv*xg);

X_v = [
   -0.5840   -0.4651
   -0.1807   -0.7064
    1.0798   -0.7788
    1.3739   -0.5938
    1.3487   -0.4330
   -0.4412    0.8298
    0.6261    0.2426
   -0.8782    0.7172
   -0.7941   -0.1354
   ];