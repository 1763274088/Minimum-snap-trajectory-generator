function xT = findTraj(r, n, m, dim, tDes, posDes)

    H_joint = [];
    
    for i=1 : m
        % find cost matrix for each segment
        H = findCostMatrix(n, r); 
        % multiply by time factor to nondimensionalize
        H = 1 ./ ((tDes(i+1, 1) - tDes(i, 1))^(2*r-1)) .* H;
        % put in block diagonal matrix
        H_joint = blkdiag(H_joint, H);
    end

    % construct equality constraints 
    [A_fixed, b_fixed] = findFixedConstraints(r, n, m, dim, posDes, tDes);
    [A_cont, b_cont] = findContConstraints(r, n, m, dim, posDes, tDes);

    % put each A_eq for each dimension into block diagonal matrix
    A_eq = [A_fixed; A_cont];
    b_eq = [b_fixed; b_cont];

    % find optimal trajectory through quadratic programming
    [xT_all, ~] = quadprog(H_joint, [], [], [], A_eq, b_eq);
    
    % explicitly break tracjetory into its piecewise parts for output
    xT = zeros((n+1), m);
    for j=1 : m
        xT(:, j) = xT_all((j-1)*(n+1)+1 : j*(n+1));
    end

end
