function [constKer,M,samples,sampleserr,error_LFR0] = P_constKer_for_LFR(matlfr,varargin)
%%   
%  File: P_constKer_for_LFR.m
%  Directory: 1_PhD_projects/00_my_toolboxes/FinslerTools/v11
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. March 10.
%  Major review on 2020. March 11. (2019b, v11)
%  Minor review on 2020. April 13. (2019b, v11)
%  Major review on 2020. April 18. (2019b, v12)
%

%%

if isa(matlfr,'lfr')
    matlfr = plfr(matlfr);
end
np = matlfr.np;

[~,~,~,~,e,f] = size(matlfr.lfrtbx_obj);
if e+f == 0
    M = matlfr.M;
    constKer = null(M);
    samples = [];
    sampleserr = [];
    error_LFR0 = norm(M * constKer);
    return
end

args.proj = [];    % Projection transformation applied to the sample points
args.tol = 1e-7;
args.N = 0;        % Nr. of samples in one iteration
args.samples = []; % Samples given preliminarily for kernel computation
args.sampleserr = []; % Samples given preliminarly for error checking
args.lims = np;       % Nr. of parameters or parameter bounds
args.extractif = false;
args = parsepropval(args,varargin{:});

DO_PROJECTION = true;
if isempty(args.proj)
    DO_PROJECTION = false;
    args.proj = @(W) W;
end

%%

TMP_RNSQSkSbcBzsqumHhrSV = pcz_dispFunctionName;

[ny,nu] = size(matlfr);

pcz_dispFunction2('Constant kernel computation for G: %dx%d, with tolerance: %g.',...
    ny,nu,args.tol)
pcz_dispFunctionSeparator;

% 2020.04.13. (április 13, hétfő), 14:58
Optimal_N = ceil(nu / ny);
N = max([args.N,Optimal_N]);
while ~isempty(find([size(matlfr) np] == N,1))
    N = N + 1;
end

M = zeros(0,nu);

if ~isempty(args.samples) && ~isempty(args.sampleserr)
    %% If the samples are given preliminarily

    warning 'NOT TESTED: 2020.04.10. (április 10, péntek), 20:50'
    
    .......................................................................
    ... MAIN COMPUTATIONS .................................................
    .......................................................................

    % Load samples for kernel computation if it is given.
    samples = args.samples;
    Nr_Samples = size(samples,2);

    M = matlfr.val(samples');
    constKer = null(M);

    pcz_dispFunction2('Kernel computed for M: %dx%d, where (%d = %d*%d)',...
        Nr_Samples*ny,nu,Nr_Samples*ny,Nr_Samples,ny);

    .......................................................................
    ... ERROR COMPUTATIONS ................................................
    .......................................................................

    % Load samples for error computation if it is given.
    sampleserr = args.sampleserr;
    N = size(sampleserr,2);

    % 2019.09.05. (szeptember  5, csütörtök), 16:17
    M = [ M ; matlfr.val(sampleserr')];
    Nr_Samples = Nr_Samples + N;

    error_LFR0 = norm(M * constKer);

    pcz_dispFunction2('Error computed for M: %dx%d, where (%d = %d*%d). Error: %d',...
        Nr_Samples*ny,nu,Nr_Samples*ny,Nr_Samples,ny, error_LFR0);

    
    .......................................................................
    ... CHECK SOLUTION ....................................................
    .......................................................................


    % norm_LFR0 = distlfr(LFR0, minlfr(LFR0*0));
    if error_LFR0 < 1e-7
        return;
    end

    pcz_dispFunction2('Constant kernel with N = %d not found, trying again with N = %d.', ...
        Nr_Samples-N, Nr_Samples+N)
    pcz_dispFunction2('Though the samples were given by the user.')

end

sampleserr = [];
error_LFR0 = [];
Nr_Samples = N;
[M,samples] = dbsample(matlfr,N,'cellout',true,'extractif',args.extractif,'proj',args.proj);
if ~isempty(M)
    M = vertcat(M{:});
end
constKer = null(M);

Nr_iterations = ceil(10 * Optimal_N / N);
i = 1;
while true
    %% If the samples are not given, nor the limits: generate random points

    pcz_dispFunction2('Kernel computed for M: (%d*%d=)%dx%d',...
        Nr_Samples,ny,size(M));

    if isempty(constKer)
        pcz_info(1,'Kernel = {0}.')
        break
    end
    .......................................................................
    ... ERROR COMPUTATION .................................................
    .......................................................................

    ZERO_lfr = minlfr(matlfr.lfrtbx_obj*constKer);
    
    [~,~,~,~,m1,~] = size(ZERO_lfr);
    if m1 == 0
        ZERO = ZERO_lfr.d;
    else
        [ZERO_sampled,~] = dbsample(plfr(ZERO_lfr,matlfr.subsvars),N,'cellout',true,'extractif',args.extractif,'proj',args.proj);
        ZERO = vertcat(ZERO_sampled{:});        
    end

    % Taking the maximum value
    error_LFR0 = max(abs(ZERO(:)));

    pcz_dispFunction2('Error computed for Merr: (%d*%d=)%dx%d. Maximum error: %d',...
        N,ny,N*ny,nu,error_LFR0);

    
    .......................................................................
    ... DISPLAY STATUS ....................................................
    .......................................................................

    pcz_dispFunctionSeparator;

    if error_LFR0 < args.tol
        pcz_info(true,'Kernel found with the given tolerance.')
        break;
    elseif i >= Nr_iterations
        pcz_info(false,'Maximum nr. iterations (%d) reached. Kernel NOT found with the given tolerance.', ...
            Nr_iterations)
        break;
    else
        i = i + 1;
    end
    
    .......................................................................
    ... MAIN COMPUTATIONS .................................................
    .......................................................................


    pcz_dispFunction2('Constant kernel with N = %d not found, trying again with N = %d.', ...
        Nr_Samples-N, Nr_Samples+N)

    [M_sampled,new_samples] = dbsample(matlfr,N,'cellout',true,'extractif',args.extractif,'proj',args.proj);
    samples = [samples new_samples];

    M = [ M ; vertcat(M_sampled{:}) ];
    Nr_Samples = Nr_Samples + numel(new_samples);

    constKer = null(M);
    
end

pcz_dispFunctionEnd(TMP_RNSQSkSbcBzsqumHhrSV);

end

function test1_simple_example
%%
    G_reset

    x = sym('x',[2 1]);

    lfrs x1 x2 
    PI = [
        1 0 
        0 1 
        x1 0
        0 x2
        x1^2*x2 0
        0 x1*x2
        x1*x2^2 x1*x2
        ];

    [m,n] = size(PI);

    % N_XI = kronLFR([1;x1;x2],eye(m));
    N_Xi = [
        eye(m)
        x1*eye(m)
        x2*eye(m)
        ];

    Aminek_keresem = N_Xi * PI;

    [N0,~,M] = P_constKer_for_LFR(Aminek_keresem','N',1);

    [U,Sigma,V] = svd(M);

    Size_N0 = size(N0)
    Size_M = size(M)

    format long g
    diag(Sigma)
    format short
    
    ZERO = plfr(N0' * Aminek_keresem);
    
    pcz_fhzero_report(ZERO,x)

end

function test2_trivial_kernel
%%

    x = sym('x',[3 1]);

    lfrs x1 x2 x3
    PI = [
        1 0 0
        0 1 0
        x1 0 1
        0 x2 x1*x3
        x1^2*x2 0 x3
        0 x1*x2 x2*x3^2
        x1*x2^2 x1*x2 1/(x1^2 + x3^2 + 1)
        ];

    [m,n] = size(PI);

    % N_XI = kronLFR([1;x1;x2],eye(m));
    N_Xi = [
        eye(m)
        x1*eye(m)
        x2*eye(m)
        ];

    Aminek_keresem = N_Xi * PI;

    [N0,~,M] = P_constKer_for_LFR(Aminek_keresem');

    [U,Sigma,V] = svd(M);

    Size_N0 = size(N0)
    Size_M = size(M)

    format long g
    diag(Sigma)
    format short

    ZERO = plfr(N0' * Aminek_keresem);
    
    pcz_fhzero_report(ZERO,x)
    
end

function test3
%%

    G_reset

    x = sym('x',[8 1]);

    lfrs x1 x2 x3 x4 x5 x6 x7 x8
    PI = [
        1 0 0
        0 x8 0
        x5*x3^2/(1 + 72^2*x1^2 + x8*x1) 0 x1*x2 
        1/(x1^2 + x3^2 + 1) x1*x2 1
        x1 0 1
        0 x2 x1*x3
        x1^2*x2*(1+x2^2*x3)/(6+x7) 0 x3
        0 x1*x5 x2*x3^2/(1 + x4^2*x1^2 + x2*x1)
        x1*x2^2 x1*x2 1/(x1^2 + x3^2 + 1)
        ];
    [m,n] = size(PI);
    
    S = randn(m+4,m);
    
    PI = S*PI;
    [m,n] = size(PI);

    % N_XI = kronLFR([1;x1;x2],eye(m));
    N_Xi = [
        eye(m)
        x1*eye(m)
        x2*eye(m)
        x3*eye(m)
        x4*eye(m)
        x5*eye(m)
        x6*eye(m)
        x7*eye(m)
        x8*eye(m)
        x1^2*eye(m)
        x1*x2*eye(m)
        x2^2*eye(m)
        x1*x3*eye(m)
        x2*x3*eye(m)
        x3^2*eye(m)
        x5^2*eye(m)
        x5*x6*eye(m)
        x6^2*eye(m)
        x5*x7*eye(m)
        x6*x7*eye(m)
        x7^2*eye(m)
        ];

    Aminek_keresem = N_Xi * PI;

    [N0_v1,~,M] = P_constKer_for_LFR(Aminek_keresem');
    
    [N0_v2,~,M] = P_constKer_for_LFR(Aminek_keresem','N',20);

    [U,Sigma,V] = svd(M);

    Size_M = size(M)

    format long g
    diag(Sigma)
    format short

    ZERO = plfr(N0_v1' * Aminek_keresem);
    pcz_fhzero_report(ZERO,x,'Elso')

    ZERO = plfr(N0_v2' * Aminek_keresem);
    pcz_fhzero_report(ZERO,x,'Masodik')

    Size_N0_v1 = size(N0_v1)
    Size_N0_v2 = size(N0_v2)

end
