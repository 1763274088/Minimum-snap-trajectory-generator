function [xT] = findTrajCorr(r, n, m, d, tDes, posDes, ineqConst)

    % construct cost matrix H
    H_opt = [];

    for dim=1 : d
        H_joint = [];
        for i=1 : m
            H = findCostMatrix(n, r); % find cost matrix for each segment
            H = 1./((tDes(i+1, 1)-tDes(i, 1))^(2*r-1)) .* H;
        
            H_joint = blkdiag(H_joint, H);
        end     
        H_opt = blkdiag(H_opt, H_joint);
    end

    % construct equality constraints
    A_opt = [];
    b_opt = [];

    for dim=1 : d
    
        % construct fixed value constraints
        [A_fixed, b_fixed] = findFixedConstraints(r, n, m, dim, posDes, ...
            tDes);
        [A_cont, b_cont] = findContConstraints(r, n, m, dim, posDes, tDes);
    
        % put each A_eq for each dimension into block diagonal matrix
        A_opt = blkdiag(A_opt, [A_fixed; A_cont]);
        b_opt = [b_opt; [b_fixed; b_cont]];   
    end

    % construct any inequality constraints
    [A, b] = constructCorrConstraints(n, m, d, posDes, ineqConst);

    % find optimal trajectory through quadratic programming
    xT_all = quadprog(H_opt,[], A, b, A_opt, b_opt);

    % explicitly break trajectory into its piecewise parts and dimensions 
    % for output
    xT = zeros((n+1), m, d);
    for dim=1 : d
        thisxT = xT_all((dim-1) * (n+1) * m + 1 : dim * (n+1) * m);
        for j=1 : m
            xT(:, j, dim) = thisxT((j-1) * (n+1) + 1 : j * (n+1));
        end
    end

end







