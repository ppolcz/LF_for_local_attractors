function X = copy(X, Y)
%% copy
%  
%  File: copy.m
%  Directory: 1_PhD_projects/00_my_toolboxes/FinslerTools/v11/@PAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. October 22.
%

%%

X.name = Y.name;
X.caption = Y.caption;
X.type = Y.type;

X = X.set_subsvars(Y.subsvars);

if Y.issym
    X = X.generate_symbolic;
end

end