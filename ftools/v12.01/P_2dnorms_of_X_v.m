function [v, ak, fci, proj_x, proj_xp, ek, Vk] = P_2dnorms_of_X_v(v)
%%
%  File: P_2dnorms_of_X_v.m
%  Directory: 7_ftools/ftools/v12.01
%  Author: Peter Polcz (ppolcz@gmail.com) 
% 
%  Created on 2021. June 02. (2021a)


%% Section 1. Imported from:
% 
%  file:   P_corner_points2.m
%  author: Polcz Péter <ppolcz@gmail.com> 
%  
%  Created on 2016.04.19. Tuesday, 11:10:42
%

I = convhull(v);
v = v(I(1:end-1), :);
Nr = size(v,1);
nx = 2;
fci = [ 1:Nr ; circshift(1:Nr,1) ]';

%% Section 2. Imported from:
% 
%  file:   P_ndnorms_of_X.m
%  author: Polcz Péter <ppolcz@gmail.com> 
%  
%  Created on 2016.04.17. Sunday, 16:40:39
%  Modified on 2019. November 19. (2019a) -- projection extended xp space
%  Modified on 2020. January 05. (2019b) -- additional segments (multiple)
%

ek = zeros(Nr,nx);
ak = zeros(Nr,nx);
Vk = cell(Nr,1);
proj_x = cell(Nr,1);
proj_xp = cell(Nr,1);

Nr_vertices_per_facet = 2^(nx-1);

for i = 1:Nr

    r0 = v(fci(i,1), :)';

    % 2D: W = [ B - A ] - actually calculated
    % 3D: W = [ B - A ; C - A ; D - A ] - actually calculated
    % 3D: W = [ B - A ; D - A ] - it would be enough
    % 4D: W = [ A2 - A1 ; B1 - A1 ; B2 - A1 ; C1 - A1 ; C2 - A1 ; D1 - A1 ; D2 - A1 ] - actually calculated
    % 4D: W = [ A2 - A1 ; B1 - A1 ; C1 - A1 ; D1 - A1 ] - it would be enough
    V = (  v(fci(i,2:Nr_vertices_per_facet),:) - repmat(r0', [Nr_vertices_per_facet-1, 1])  )';
    
    % normal vector (norm(n) == 1)
    nvec = null(V');

    % [3.1] ek
    ek(i,:) = nvec;

    % distance of the origin to the dim dimensional facet
    dst = dot(nvec, r0);
    
    % [3.2] ak (may be undefined - division by 0)
    ak(i,:) = nvec / dst;
    
    % [3.3] Vk
    V = orth(V);
    Vk{i} = V;
    
    % [3.4.1] projection (only x needed to be substituted)
    proj_x{i} = @(x) V*V'*x + dst*nvec*ones(1,size(x,2));
    
    % V_xp = blkdiag(V,eye(np));
    % nvec_xp = [nvec ; zeros(np,1)];
    % 
    % % [3.4.2] projection in extended space ([x;p] needed to be substituted)
    % proj_xp{i} = @(xp) V_xp*V_xp'*xp + dst*nvec_xp*ones(1,size(xp,2));
end


end


function test1
    %%
    
    X_v = [
       -0.2521    1.0387
       -1.4218    0.0202
       -0.5849   -1.2706
        1.0487   -0.8471
        1.3008    0.8370
        ];
    
    [v,ak,fci,proj] = P_2dnorms_of_X_v(X_v)
    
    
    for i = 1:size(fci,1)
    
        W = randn(2,50);
        Wi = proj{i}(W);
    
        pcz_plotpoint(Wi), hold on;
        
    end
    
    axis equal
    axis tight
    
end

function test2
    %%
    
    x_lim = [
        -1 1 
        -1 1
        -1 1
        ];
    
    [v,ak,fci,proj] = P_ndnorms_of_X(x_lim)
    
    
    for i = 1:size(fci,1)
    
        W = randn(3,50);
        Wi = proj{i}(W);
    
        pcz_plotpoint(Wi), hold on;
        
    end
    
    axis equal
    axis tight
    
end