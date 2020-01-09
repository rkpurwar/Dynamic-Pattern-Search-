% Computes motion vectors using Dynamic Pattern Search method based on spatial and temporal motion vector correlations
% Input
%   imgP : The image for which we want to find motion vectors
%   imgI : The reference image
%   mbSize : Size of the macroblock
%   p : Search parameter  (read literature to find what this means)
%   temporal_adj_mv : temporal motion vector 
%
% Ouput
%   motionVect : the motion vectors for each integral macroblock in imgP
%   DPScomputations: The average number of points searched for a macroblock



function [motionVect, DPScomputations] = motionEstDPS(imgP, imgI, mbSize, p, temporal_adj_mv)

[row col] = size(imgI);
tic
vectors = zeros(2,row*col/mbSize^2);
costs = ones(1, 7) * 65537;


% The index points for Small Diamond search pattern
SDSP(1,:) = [ 0 -1];
SDSP(2,:) = [-1  0];
SDSP(3,:) = [ 0  0];
SDSP(4,:) = [ 1  0];
SDSP(5,:) = [ 0  1];

% We will be storing the positions of points where the checking has been
% already done in a matrix named checkMatrix[]. This matrix is initialised to zero and as a point is
% checked, we set the corresponding element in the matrix to one. 

checkMatrix = zeros(2*p+1,2*p+1);

computations = 0;


% we start off from the top left of the image
% we will walk in steps of mbSize
% mbCount will keep track of how many blocks we have evaluated

mbCount = 1;
for i = 1 : mbSize : row-mbSize+1
    for j = 1 : mbSize : col-mbSize+1
        
          
        x = j;
        y = i;
     
        costs(3) = costFuncMAD(imgP(i:i+mbSize-1,j:j+mbSize-1), ...
                                    imgI(i:i+mbSize-1,j:j+mbSize-1),mbSize);
        
        checkMatrix(p+1,p+1) = 1;
        computations =  computations + 1; 
        % if we are in the left most column then we have to make sure that
        % we just do the Large Diamond Search Pattern (LDSP) with stepSize = 2
        if (j-1 < 1)
            stepSize = 2;
            maxIndex = 5;
        else 
            stepSize1 = max(abs(vectors(1,mbCount-1)), abs(vectors(2,mbCount-1)));
            stepSize2 = max(abs(temporal_adj_mv(1,mbCount)),abs(temporal_adj_mv(2,mbCount)));
            stepSize = max(stepSize1,stepSize2);

            % now we have to make sure that if the point due to motion
            % vector is one of the LDSP points then we dont calculate it
            % again
            if ( (abs(vectors(1,mbCount-1)) == stepSize && vectors(2,mbCount-1) == 0) ...
                 || (abs(vectors(2,mbCount-1)) == stepSize && vectors(1,mbCount-1) == 0) ...
                 || (abs(temporal_adj_mv(1,mbCount)) == stepSize && temporal_adj_mv(2,mbCount) == 0) ...
                 || (abs(temporal_adj_mv(2,mbCount)) == stepSize && temporal_adj_mv(1,mbCount) == 0) ) ...
                 
                maxIndex = 5; % we just have to check at the rood pattern 5 points
                
            else
                maxIndex = 7; % we have to check 6 points
                LDSP(6,:) = [ vectors(2, mbCount-1)  vectors(1, mbCount-1)];
                LDSP(7,:) = [ temporal_adj_mv(2, mbCount)  temporal_adj_mv(1, mbCount)];
            end
        end
        
        % The index points for first and only Large Diamond search pattern
        
        LDSP(1,:) = [ 0 -stepSize];
        LDSP(2,:) = [-stepSize 0]; 
        LDSP(3,:) = [ 0  0];
        LDSP(4,:) = [stepSize  0];
        LDSP(5,:) = [ 0  stepSize];
        
        
        % do the LDSP
        
        
        for k = 1:maxIndex
            refBlkVer = y + LDSP(k,2);   % row/Vert co-ordinate for ref block
            refBlkHor = x + LDSP(k,1);   % col/Horizontal co-ordinate
            if ( refBlkVer < 1 || refBlkVer+mbSize-1 > row ...
                 || refBlkHor < 1 || refBlkHor+mbSize-1 > col)
             
                continue; % outside image boundary
            end

            if (k == 3 || stepSize == 0)
                continue; % center point already calculated
            end

            costs(k) = costFuncMAD(imgP(i:i+mbSize-1,j:j+mbSize-1), ...
                  imgI(refBlkVer:refBlkVer+mbSize-1, refBlkHor:refBlkHor+mbSize-1), mbSize);
            computations =  computations + 1;
            checkMatrix(LDSP(k,2) + p+1, LDSP(k,1) + p+1) = 1;
            
        end
        
        [cost, point] = min(costs);
        
        
        % The doneFlag is set to 1 when the minimum
        % is at the center of the diamond           

        x = x + LDSP(point, 1);
        y = y + LDSP(point, 2);
        costs = ones(1,5) * 65537;
        costs(3) = cost;
       

        doneFlag = 0;   
        if ((point == 6) || (point == 7))
            doneFlag = 1;
        end 
        while (doneFlag == 0)
            for k = 1:5
                refBlkVer = y + SDSP(k,2);   % row/Vert co-ordinate for ref block
                refBlkHor = x + SDSP(k,1);   % col/Horizontal co-ordinate
                if ( refBlkVer < 1 || refBlkVer+mbSize-1 > row ...
                      || refBlkHor < 1 || refBlkHor+mbSize-1 > col)
                      continue;
                end

                if (k == 3)
                    continue
                elseif (refBlkHor < j-p || refBlkHor > j+p || refBlkVer < i-p ...
                            || refBlkVer > i+p)
                        continue;
                elseif (checkMatrix(y-i+SDSP(k,2)+p+1 , x-j+SDSP(k,1)+p+1) == 1)
                    continue
                end

                costs(k) = costFuncMAD(imgP(i:i+mbSize-1,j:j+mbSize-1), ...
                             imgI(refBlkVer:refBlkVer+mbSize-1, ...
                                 refBlkHor:refBlkHor+mbSize-1), mbSize);
                checkMatrix(y-i+SDSP(k,2)+p+1, x-j+SDSP(k,1)+p+1) = 1;
                computations =  computations + 1;
                
  
            end
            
            [cost, point] = min(costs);
           
            if (point == 3)
                doneFlag = 1;
            else
                x = x + SDSP(point, 1);
                y = y + SDSP(point, 2);
                costs = ones(1,5) * 65537;
                costs(3) = cost;
            end
            
        end  % while loop ends here
        
        vectors(1,mbCount) = y - i;    % row co-ordinate for the vector
        vectors(2,mbCount) = x - j;    % col co-ordinate for the vector            
        mbCount = mbCount + 1;
        costs = ones(1,7) * 65537;
        
        checkMatrix = zeros(2*p+1,2*p+1);
    end
end
    
motionVect = vectors;
DPScomputations = computations/(mbCount-1) ; 
 toc   
    
 