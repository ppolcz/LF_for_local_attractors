function [A,B,C,D,PI,PI_1] = gss2pgen(sysgss)
%%
%  File: lfr2pgen.m
%  Directory: 7_ftools/ftools/v12.01/utilities/gss_utils
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2021. January 20. (2020b)
%

if ~isa(sysgss,'gss')
    error 'We require a `gss` object!'
end

% CALL
% [a,bw,bu,cz,cy,dzw,dzu,dyw,dyu,st,no,ni,ndo,ndi,ns] = gssdata(sys)
%
% OUTPUT ARGUMENTS
% - a,bw,bu,cz,cy,dzw,dzu,dyw,dyu: state-space matrices of sys.
% - st: sample time  in seconds (-1 if unspecified, 0 if M is a con-
%       tinuous-time model).
% - no: number of external outputs (size of y).
% - ni: number of external inputs (size of u).
% - ndo: number of Delta-block outputs (size of w).
% - ndi: number of Delta-block inputs (size of z).
% - ns: number of states of M.
[~,~,~,~,~,D,C,B,A,~,~,nu,m1,~,~] = gssdata(sysgss);
Delta = sysgss.D;

% sysgss.M = ss([],[],[],[ D C ; B A ])
% sysgss.D = "Delta" (array of structs corresponding to each parameter)

m = nu + m1;

M = ss([],[],[],[D C ; eye(m)]);
PI = gss(M,Delta);

M_1 = ss([],[],[],[D C ; eye(m1,m)]);
PI_1 = gss(M_1,Delta);

end