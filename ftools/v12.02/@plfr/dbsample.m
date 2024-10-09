function [sampled,samples_subsvars,samples_lfrvars,samples_good] = dbsample(syslfr,varargin)
%%
%  File: dbsample.m
%  Directory: 7_ftools/ftools/v12/@plfr
%  Author: Peter Polcz (ppolcz@gmail.com)
%
%  Created on 2020. May 19. (2019b)
%
% Instructions for plfr/dbsample:
%
%  N:    number of samples (default 1)
%
%  proj: project samples (default @(x) x)
%        e.g. @(x) V*V'*x
%        where x = [ sample1.x1 sample2.x1 ... sampleN.x1
%                    ...        ...            ...
%                    sample1.xn sample2.xn ... sampleN.xn ]
%        [x1;x2;...;xn] corresponds to subsvars of syslfr (by default)
%
%  lims: bounds for x (default: computed from syslfr)
%        e.g. lims = [ x1.min x1.max
%                      ...    ...
%                      xn.min xn.max ];
%        [x1;x2;...;xn] corresponds to subsvars of syslfr (by default)
%
%  lfrspace: proj and lims corresponds to lfrvars  (if true)
%            proj and lims corresponds to subsvars (if false: default)

%% SYSLFR is converted to plfr

if isa(syslfr,'gss')
    syslfr = gss2lfr(syslfr);
end

if isa(syslfr,'lfr')
    syslfr = plfr(syslfr);
end

%% Pares arguments
% 
% args.lims = [];
% args.proj = [];
% args.lfrspace = false;
% args.N = 3;
% args.verbose = false;
% args.cellout = false;
% args = parsepropval(args,varargin{:});



p = inputParser;
addOptional(p,'N',1,@(x) isnumeric(x) && isscalar(x) && (x >= 1));
addParameter(p,'proj',[]);
addParameter(p,'lims',[]);
addParameter(p,'lfrspace',false);
addParameter(p,'verbose',false);
addParameter(p,'cellout',false);
addParameter(p,'extractif',false);

parse(p,varargin{:});
args = p.Results;

%% LIMS corresponds to the lfr variables (not the substitutional variables)

if isempty(args.lims)
    lims = syslfr.bounds;
elseif args.lfrspace
    lims = args.lims;
else
    lims = syslfr.subsvars2lfrvars_fh(args.lims);
end

MAX = max(lims(~isinf(lims))) * 10;
MIN = min(lims(~isinf(lims))) * 10;
lims(isinf(lims) & lims < 0) = MIN;
lims(isinf(lims)) = MAX;

%% PROJ corresponds to the lfr variables (not the substitutional variables)

if isempty(args.proj)
    proj = @(x) x;
elseif args.lfrspace
    proj = args.proj;
else
    proj = @(x) syslfr.subsvars2lfrvars_fh(args.proj(syslfr.lfrvars2subsvars_fh(x)));
end

N = args.N;

%%

samples_cell = cellfun(@(a,b) {rand(1,N)*(b-a)+a}, num2cell(lims(:,1)), num2cell(lims(:,2)));
samples_raw = vertcat(samples_cell{:});
samples_lfrvars = proj(samples_raw);
samples_good = prod(lims(:,1)*ones(1,N) <= samples_lfrvars & samples_lfrvars <= lims(:,2)*ones(1,N),1) == 1;

if args.extractif
    samples_lfrvars = samples_lfrvars(:,samples_good);
    samples_good = samples_good(:,samples_good);
end

samples_subsvars = syslfr.lfrvars2subsvars_fh(samples_lfrvars);

%%

if args.cellout
    sampled = syslfr.val('cellout',samples_subsvars);
else
    sampled = syslfr.val(samples_subsvars);
end

%%

if args.verbose
    % 1: reset SCOPE_DEPTH if not set
    % 2: do not modify verbosity (VERBOSE) if not set (it is TRUE by default)
    G_reset(12);

    TMP_zEmxYCBmzGJutZbUVPfn = pcz_dispFunctionName;

    stringify = @(msg,s) [ msg ':' newline '[ ' cell2mat(join(cellfun(@char, num2cell(s), 'UniformOutput', 0),[newline '  '])) ' ]' newline ];

    subsvars = sym('w',[syslfr.np 1]);
    lfrvars = syslfr.subsvars2lfrvars_fh(subsvars);
    projected = proj(lfrvars);

    pcz_dispFunction2('Variable bounds in the LFR''s variable space:')
    pcz_dispFunction_num2str(lims)
    pcz_dispFunction2('\nNumber of random samples: %d\n', N)
    pcz_dispFunction2(stringify('Substitutional variables',subsvars))
    pcz_dispFunction2(stringify('LFR variables', lfrvars))
    pcz_dispFunction2(stringify('Projected vector',projected))

    pcz_dispFunction_num2str(samples_raw)
    pcz_dispFunction_num2str(samples_lfrvars)
    pcz_dispFunction_num2str(samples_subsvars)

    if args.cellout
        First_sample = sampled{1};
    else
        First_sample = sampled(1:syslfr.ny,1:syslfr.nu);
    end 
    pcz_dispFunction_num2str(First_sample);

    pcz_dispFunctionEnd(TMP_zEmxYCBmzGJutZbUVPfn);
end


end
