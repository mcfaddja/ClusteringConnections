function [ newA ] = SVDgaps( A, k, testRow, testCol, terms, objcs )

%% INPUT: A = m by n matrix.
%% k = approximation order for SVD.
%% testRow = value to test the row scores against.
%% testCol = value to test the column scores against.
%% terms = row labels
%% objcs = column labels

% Size info for A.
sizeA = size(A);
m = sizeA(1,1);
n = sizeA(1,2);





% Approximate SVD of order k for A.
[Uk, Sk, Vk] = svds(A, k); 


% Clustering rows.
[sortedU, UkIndex] = sort(Uk); % sort Uk vectors
UkGapMatrix = sortedU(2:m, :) - sortedU(1:m-1, :); % calculate gaps between adjacent Uk elements
UkGapMeans = mean(UkGapMatrix, 1); % find mean gap size for Uk elements
UkGapStDev = std(UkGapMatrix, 1, 1); % find std of gap size for Uk elements
UkGapScore = (UkGapMatrix - ones(m - 1, 1) * UkGapMeans) ./ (ones( m - 1, 1) * UkGapStDev); % score the Uk gaps
[UkRow, UkCol] = find(UkGapScore > testRow); 
Uk0 = sparse( UkRow, UkCol, ones( length(UkRow), 1 ), m, k);

UkC = sparse(m, k);
for j = 1:k
    cnt = 0;
    
    for i = 1:m
        UkC(i,j) = cnt + 1; % creates cluster labels
        if full( Uk0(i,j) ) == 1
            cnt = cnt + 1; % cluster label changes 
        end
    end
end

[SortedIndexU, RowIndex] = sort(UkIndex, 1);
for i = 1:k
    UkC(:,i) = UkC( (RowIndex(:,i)), i);
end

[UkB, UkI, UkH] = unique(UkC, 'rows'); % finds rows with the same label patterns
[termsClust, termsClustInd] = sort(UkH);
UkC(termsClustInd, :);
termsClustNum = size(UkB);



% Clustering cols.
[sortedV, VkIndex] = sort(Vk); % sort Vk vectors
VkGapMatrix = sortedV(2:n, :) - sortedU(1:n-1, :); % calculate gaps between adjacent Vk elements
VkGapMeans = mean(VkGapMatrix, 1); % find mean gap size for Vk elements
VkGapStDev = std(VkGapMatrix, 1, 1); % find std of gap size for Vk elements
VkGapScore = (VkGapMatrix - ones(n - 1, 1) * VkGapMeans) ./ (ones( n - 1, 1) * VkGapStDev); % score the Vk gaps
[VkRow, VkCol] = find(VkGapScore > testCol); 
Vk0 = sparse( VkRow, VkCol, ones( length(VkRow), 1 ), n, k);

VkD = sparse(n,k);
for j = 1:k
    cnt = 0;
    
    for i = 1:n
        VkD(i,j) = cnt + 1; % creates cluster labels
        if full( Vk0(i,j) ) == 1
            cnt = cnt + 1; % cluster label changes 
        end
    end
end

[SortedIndexV, ColIndex] = sort(VkIndex, 1);
for i = 1:k
    VkD(:,i) = VkD( (ColIndex(:,i)), i);
end

[VkB, VkI, VkH] = unique(VkD, 'rows'); % finds cols with the same label patterns
[objcsClust, objcsClustInd] = sort(VkH);
VkD(objcsClustInd, :);
objcsClustNum = size(VkB);



newA = A(termsClustInd, objcsClustInd);
newTerms = terms(termsClustInd);
newObjcs = objcs(objcsClustInd);
















end

