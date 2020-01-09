% Computes motion vectors using Modified Dynamic Pattern Search method based on spatial and temporal motion vector correlations
% Input
%   imgP : The image for which we want to find motion vectors
%   imgI : The reference image
%   mbSize : Size of the macroblock
%   p : Search parameter  (read literature to find what this means)
%
% Ouput
%   motionVect : the motion vectors for each integral macroblock in imgP
%   MDPScomputations: The average number of points searched for a macroblock


function [motionVect, MDPScomputations] = motionEstMDPS(imgP, imgI, mbSize, p, temporal_adj_mv)
[row col] = size(imgI);
vectors = zeros(2,row*col/mbSize^2);
costs = ones(1, 5) * 65537;
adjacent_cost = ones(1,2)*65537;

% The index points for Small Diamond search pattern
SDSP(1,:) = [ 0 -1];
SDSP(2,:) = [-1  0];
SDSP(3,:) = [ 0  0];
SDSP(4,:) = [ 1  0];
SDSP(5,:) = [ 0  1];

% We will be storing the positions of points where the checking has been
% already done in a matrix named CheckMatrix[]. This matrix is initialised to zero and as a point is
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
         
         
         if (j-1 >= 1) 
             temporal_mag = sqrt(sum(temporal_adj_mv(:,mbCount).*temporal_adj_mv(:,mbCount)));  % find magnitude of temporally coherent vector
             spatial_mag = sqrt(sum(vectors(:,mbCount-1).*vectors(:,mbCount-1)));               % find magnitude of spatially coherent vector
         
                 if (temporal_mag ~= 0) && (spatial_mag ~= 0)                                   % when both adjacent vectors are non zero
                   
                 % find angle of cosine between spatial and temporal
                 % vectors
                    angle = acos(sum(temporal_adj_mv(:,mbCount).* vectors(:,mbCount-1))/(temporal_mag*spatial_mag))*180/pi; 
                         
                     refBlkVer1 = y + temporal_adj_mv(1,mbCount);                               % row/Vert co-ordinate for ref block using tempral vector
                     refBlkHor1 = x + temporal_adj_mv(2,mbCount);                               % col/Horizontal co-ordinate for ref block using temporal vector
                     refBlkVer2 = y + vectors(1,mbCount-1);                                     % row/Vert co-ordinate for ref block using spatial vector
                     refBlkHor2 = x + vectors(2,mbCount-1);                                     % col/Horizontal co-ordinate for ref block using spatial vector
                    if ( refBlkVer1 < 1 || refBlkVer2 < 1 || refBlkVer1+mbSize-1 > row || refBlkVer2+mbSize-1 > row ...
                    || refBlkHor1 < 1 || refBlkHor2 < 1 || refBlkHor1+mbSize-1 > col || refBlkHor2+mbSize-1 > col)
             
                        continue; % outside image boundary
                    
                    end 
        
         
         if (angle <= 45)   % if the angle between spatial & temporal vectors is less than pi/4  
             
                        % compute cost at the location of temporal vector
                        adjacent_cost(1) = costFuncMAD(imgP(i:i+mbSize-1,j:j+mbSize-1), ...
                        imgI(refBlkVer1:refBlkVer1+mbSize-1, refBlkHor1:refBlkHor1+mbSize-1), mbSize);
                        computations =  computations + 1;
                        
                        % compute cost at the location of spatial vector
                        adjacent_cost(2) = costFuncMAD(imgP(i:i+mbSize-1,j:j+mbSize-1), ...
                        imgI(refBlkVer2:refBlkVer2+mbSize-1, refBlkHor2:refBlkHor2+mbSize-1), mbSize);
                        computations =  computations + 1;
                        [cost, point] = min(adjacent_cost);
                        if (point == 1)
                            vectors(:,mbCount) = temporal_adj_mv(:,mbCount);
                        else 
                            vectors(:,mbCount) = vectors(:,mbCount-1);
                        end
                           adjacent_cost = ones(1,2)*65537;
                           mbCount = mbCount + 1;
                       continue; 
                end 
            end 
         end 
        
        


             costs(3) = costFuncMAD(imgP(i:i+mbSize-1,j:j+mbSize-1), ...
                                    imgI(i:i+mbSize-1,j:j+mbSize-1),mbSize);
        
        checkMatrix(p+1,p+1) = 1;
        computations =  computations + 1; 
        % if we are in the left most column then we have to make sure that
        % we just do the Large Diamond search pattern (LDSP) with stepSize = 2
        if (j-1 < 1)
        stepSize = 2;
        else 
            stepSize1 = max(abs(vectors(1,mbCount-1)), abs(vectors(2,mbCount-1)));
            stepSize2 = max(abs(temporal_adj_mv(1,mbCount)),abs(temporal_adj_mv(2,mbCount)));
            stepSize = max(stepSize1,stepSize2);
        end 
        maxIndex = 5;
        
          % The index points for first and only LDSP
        
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
        costs = ones(1,5) * 65537;
        
        checkMatrix = zeros(2*p+1,2*p+1);
         
    end
end    
motionVect = vectors;
MDPScomputations = computations/(mbCount-1) ; 
   
        
                        
                
                 
                 
        
        
        