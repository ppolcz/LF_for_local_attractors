classdef plpv

properties
    fn
    
    dp_scalar
end

methods
    
    function o = plpv(fn)
        o.fn = fn;
        
        np = 0;
        while true
            p_cell = num2cell(zeros(np,1));
            try
                o.fn(p_cell{:});
                
        end
        
    end
    
    function r = roundOff(obj)
     r = round([obj.Value],2);
    end
    
    function r = multiplyBy(obj,n)
     r = [obj.Value] * n;
    end
end
end