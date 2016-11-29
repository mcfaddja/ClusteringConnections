function [ newA ] = SVDsigns( A, k, terms, objcs )

%% INPUT: A = m by n matrix.
%% k = approximation order for SVD.
%% terms = row labels
%% objcs = column labels

% Size info for A.
sizeA = size(A);
m = sizeA(1,1);
n = sizeA(1,2);


% Approximate SVD of order k for A.
[Uk, Sk, Vk] = svds(A, k); 


% Clustering rows.
Uk0 = (Uk >= 0); % Finds positive entries of the U matrix from the SVD.
x = zeros(m, 1);
for i = 1:k
    x = x + (2^(i - 1)) * (Uk0(:, k - i + 1)); % sign pattern of each row in U
end

[sortedRowX, rowIndex] = sort(x); % sorts x by the sign pattern
rowClustCNT = length(unique(x)); 


% Clustering Columns.
Vk0 = (Vk >= 0); % Finds positive entries of the V matrix from the SVD.
y = zeros(n, 1);
for i = 1:k
    y = y + (2^(i - 1)) * (Vk0(:, k - i + 1)); % sign pattern of each column in V
end

[sortedColY, colIndex] = sort(y); % sorts y by the sign pattern
colClustCNT = length(unique(y)); 


newA = A(rowIndex, colIndex);
newTerms = terms(rowIndex);
newObjcs = objcs(colIndex);





































end

