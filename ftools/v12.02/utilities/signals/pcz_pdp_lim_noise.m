function [t,p,dp] = pcz_pdp_lim_noise(t, p_lim, dp_lim)
%% pcz_pdp_lim_noise
%  
%  File: pcz_pdp_lim_noise.m
%  Directory: 2_demonstrations/lib/matlab
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. February 20.
%

%%

% mdl = 'sim_bounded_varbounded_noise';
% load_system(mdl);
% 
% set_param(mdl, 'SimulationMode', simMode);
% 
% simout = sim(mdl, 'StopTime', num2str(T));
% 
% 
% t = simout.tout;
% p = permute(simout.yout{1}.Values.Data,[3,1,2]);
% dp = permute(simout.yout{2}.Values.Data,[3,1,2]);


dim = size(p_lim,1);

N = numel(t);

dt = diff(t);
dt = [dt(:)' dt(end)];

dp_lim1 = dp_lim(:,1)*dt;
dp_lim2 = dp_lim(:,2)*dt;



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

t = t(:)';
p = p';
dp = dp';

end