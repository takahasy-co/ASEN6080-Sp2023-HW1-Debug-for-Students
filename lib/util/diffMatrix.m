function [deltaM, deltaMPct] = diffMatrix(m1,m2)

deltaM     = zeros(size(m1));
deltaMPct  = zeros(size(m1));
numParams1 = length(m1(:,1));
numParams2 = length(m1(1,:));

for ii = 1:numParams1
    
    for jj = 1:numParams2
        
        deltaM(ii,jj) = m2(ii,jj) - m1(ii,jj);
        
        if m1(ii,jj) ~= 0
            
            deltaMPct(ii,jj) = ( m2(ii,jj) - m1(ii,jj) ) / m1(ii,jj)*100;
            
        else
            
            deltaMPct(ii,jj) = m2(ii,jj);
            
        end % For if
        
    end % For jj
    
end % For ii
