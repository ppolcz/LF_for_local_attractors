function Monomials = pcz_monomials(x,degrees)
arguments
    x (:,1) sym
    degrees (1,:)
end


Monomials = sym([]);

all_degrees = degrees(1):degrees(end);

Phi = x;
for d = all_degrees

    if d == 0
        Monomials = [
            Monomials
            1
            ];
    elseif d == 1
        Monomials = [
            Monomials
            x
            ];
    else
        
        Phi = kron(x,eye(numel(Phi)))*Phi;
        [~,Idx] = unique(Phi);
        Phi = Phi(sort(Idx));
        
        if ismember(d,degrees)
            Monomials = [
                Monomials
                Phi
                ];
        end
        
    end
    

end

end

      

function test
%%

x = sym('x',[3 1])

m = pcz_monomials(x,0:3)

m = pcz_monomials(x,[1 3])

end
