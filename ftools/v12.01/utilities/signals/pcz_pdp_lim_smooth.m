function [t,p,dp] = pcz_pdp_lim_smooth(t, p_lim, dp_lim, varargin)
%% pcz_pdp_lim_smooth
%  
%  File: pcz_pdp_lim_smooth.m
%  Directory: 2_demonstrations/lib/matlab
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. February 27.
%

%%

opts = struct;
opts.wsize = 101;
opts = parsepropval(opts,varargin{:});


tmin = t(1);
tmax = t(end);
dt = (tmax-tmin)/(numel(t)-1);

tnew = (tmin:dt:5*tmax)';


dim = size(p_lim,1);

N = numel(tnew);

dt = diff(tnew);
dt = [dt(:)' dt(end)];


scale_factor = 10;

dp_lim1 = dp_lim(:,1)*dt*scale_factor;
dp_lim2 = dp_lim(:,2)*dt*scale_factor;



dp_normalized = min(max(0.5*randn(dim,N),-1),1);
dp = 0.5*(dp_normalized + 1).*(dp_lim2-dp_lim1) + dp_lim1;

p = zeros(dim,N);
p(:,1) = rand(dim,1).*(p_lim*[-1;1]) + p_lim(:,1);

for i = 2:N
   p(:,i) = min(max(p(:,i-1) + dp(:,i-1),p_lim(:,1)),p_lim(:,2));
end

dp = diff(p,[],2);
dp = [zeros(dim,1) dp];
dp = dp ./ (ones(dim,1)*dt);

p = p';
dp = dp';

if isscalar(opts.wsize)
    opts.wsize = ones(1,dim)*opts.wsize;
end

minlen = numel(tnew);
p_cell = cell(1,dim);
for d = 1:dim
    w = hamming(opts.wsize(d));
    w = w / sum(w);
    
    p_cell{d} = conv(p(:,d),w,'valid');

    minlen = min(minlen,numel(p_cell{d}));
end

p_cell = cellfun(@(c) {c(1:minlen,:)}, p_cell);

p = [p_cell{:}];
dp = diff(p) ./ ( diff(tnew(1:size(p,1))) * ones(1,dim) );

p = p(1:end-1,:);
t = tnew(1:size(p,1),:);

p = p(t <= tmax,:);
dp = dp(t <= tmax,:);
t = t(t <= tmax,:);

end