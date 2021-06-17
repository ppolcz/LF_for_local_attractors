function save_variables(fname,varargin)
%%
%  File: save_variables.m
%  Directory: 0_VEDETT/28_2021_03_02_attractoros_cikkbe
%  Author: Peter Polcz (ppolcz@gmail.com)
%
%  Created on 2021. April 23. (2020b)
%

s = struct;

for i = 1:2:numel(varargin)

    var = varargin{i};
    msg = varargin{i+1};
    varname = inputname(i+1);
    vartype = class(var);



    switch vartype

        case 'plfr'

            A = var.A;
            B = var.B;
            C = var.C;
            D = var.D;
            I_m = eye(var.m1);
            
            r = var.blk.desc(1,:);
            v = var.vars;

            % cellfun(@(i) {assign(sprintf('I_r%d = eye(r(%d));',i,i))}, num2cell(1:var.np));

            s.(varname).A = A;
            s.(varname).B = B;
            s.(varname).C = C;
            s.(varname).D = D;
            s.(varname).Delta.r = r;
            s.(varname).Delta.vars = var.blk.names;
            s.(varname).m = var.m1;
            s.(varname).msg = msg;
            s.(varname).lfr_formula = 'A + B / (I_m - Delta*D) * Delta*C, where Delta = SUM_i ( var_i * I_ri )';


        case 'sym'
            
            s.(varname).symbolic = var;
            s.(varname).msg = msg;
            
        otherwise

            warning('Class `%s`(%s) not supported yet',vartype,varname)

    end

end

save(fname,'-append','-struct','s')

end
