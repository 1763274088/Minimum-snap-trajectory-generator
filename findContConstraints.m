function [A_eq, b_eq] = findContConstraints(r, n, m, dim, posDes, tDes)

    A_eq = [];
    b_eq = [];

    derCoeff = findDerivativeCoeff(n, r);

    for j=0 : m       % for each keyframe
        for i=0 : r-1 % for all derivatives from 0 to r-1
            if posDes(i+1, j+1, dim) == 255, %if unfixed
                % construct and add constraint
                A_temp = [];
                b_temp = [];
            
                % if constraint is at an intermediate keyframe
                if j>0 && j<m 
                    A_temp = zeros(1, (n+1)*m);
                    maxPower = nnz(derCoeff(i+1, :)) - 1;
                
                    for k=0 : maxPower
                    
                        A_temp(1, (j-1)*(n+1)+k+1) = derCoeff(i+1, k+1);
                        A_temp(1, (j-1)*(n+1)+k+1) = 1 / ((tDes(j+1, 1) ...
                            - tDes(j, 1))^i) * A_temp(1, (j-1)*(n+1)+k+1); 
                    
                        A_temp(1, j*(n+1)+k+1) = - 0^(maxPower - k) * ...
                            derCoeff(i+1, k+1);
                        A_temp(1, j*(n+1)+k+1) = 1/((tDes(j+2, 1) - ...
                            tDes(j+1, 1))^i) * A_temp(1, j*(n+1)+k+1); 
                    end
                    b_temp = 0;
                end     
                
                A_eq = [A_eq; A_temp];
                b_eq = [b_eq; b_temp];
            end
        end
    end
end