function H = findCostMatrix(n, r)

    H = zeros(n+1);
    tempTerms = ones(n+1, 1);

    for i=0 : n
        tempTerm = 1;
        for m=0 : r-1
            tempTerm = tempTerm * (i-m);
        end
        tempTerms(i+1, 1) = tempTerm;
    end

    for i=0 : n
        for j=0 : n
            if i>=r && j>=r
                H(i+1, j+1) = tempTerms(i+1, 1) * tempTerms(j+1, 1) / ...
                    (i+j-2*r+1);
            end
        end
    end

    H = rot90(rot90(H));
    
end