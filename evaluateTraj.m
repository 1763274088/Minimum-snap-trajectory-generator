function [dxT, derivativesX] = evaluateTraj(t, n, m, d, xT, tDes, r, ...
    derivativesX)

    dxT = zeros(r + 1, d);

    if isempty(derivativesX)
        derivativesX = cell(r + 1, 1);
        derivativesX{1} = xT;

        for l=1 : r
            thisDer = zeros(n-l+1, m, d);
            for j=1 : m
                for k=1 : d
            
                    [a, ~, ~] = size(thisDer); 
                    [~, b2, ~] = size(polyder(derivativesX{l}(:, j, k)));
            
                    thisDer(:, j, k) = [zeros(1, a-b2) ...
                        polyder(derivativesX{l}(:, j, k))];
                end
            end
            derivativesX{l+1} = thisDer;
        end
    end

    % evaluate trajectory at given time
    if t < tDes(1, 1)
        % evaluate in each dimension
        for k=1 : d
            % evaluate each derivative at the first trajectory's inital 
            % time
            for l=0 : r 
                dxT(l+1, k) = 1/((tDes(2, 1) - tDes(1, 1))^(l)) * ...
                    polyval(derivativesX{l+1}(:, 1, k), 0);
            end
        end
    
        % find which piece of the trajectory we're on based on time and
        % evaluate there
    elseif t < tDes(m+1, 1)
    
        for j=1 : m
            if t < tDes(j+1, 1)
                % find the nondimensionalized time
                scaledt = (t-tDes(j, 1)) / (tDes(j+1, 1)-tDes(j, 1)); 
            
                % evaluate in each dimension
                for k=1 : d
                    for l=0 : r
                        dxT(l+1, k) = 1 / ((tDes(j+1, 1) - ...
                            tDes(j, 1))^l) * ...
                            polyval(derivativesX{l+1}(:, j, k), scaledt);
                    end
                end    
                break;
            end
        end
    
        % if after the final keyframe time, assume hover at final keyframe 
        % position
    else
    
        % evaluate in each dimension
        for k=1 : d
            for l=0 : r % evaluate each derivative 
                dxT(l+1, k) = 1 / ((tDes(m+1, 1) - tDes(m, 1))^l) * ...
                    polyval(derivativesX{l+1}(:, m, k), 1);
            end
        end
    end
end


