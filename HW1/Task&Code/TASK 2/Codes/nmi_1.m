function NMI = nmi_1( A, B ) 

if length( A ) ~= length( B)
    error('length( A ) must == length( B)');%判断前后两个分类结果的数据点数应相同，若不同则直接结束
end
total = length(A);
A_ids = unique(A);
B_ids = unique(B);%取出A和B中一样的元素（无重复）

% Mutual information
MI = 0;
for idA = A_ids
    for idB = B_ids
         idAOccur = find( A == idA );
         idBOccur = find( B == idB );
         idABOccur = intersect(idAOccur,idBOccur); 
         
         px = length(idAOccur)/total;
         py = length(idBOccur)/total;
         pxy = length(idABOccur)/total;
         
         MI = MI + pxy*log2(pxy/(px*py)+eps); % eps : 修正值，最小正数；此处先计算MI值

    end
end

% Normalized Mutual information
Hx = 0; % Entropies
for idA = A_ids
    idAOccurCount = length( find( A == idA ) );
    Hx = Hx - (idAOccurCount/total) * log2(idAOccurCount/total + eps);%计算X所蕴含的自信息（熵）
end
Hy = 0; % Entropies
for idB = B_ids
    idBOccurCount = length( find( B == idB ) );
    Hy = Hy - (idBOccurCount/total) * log2(idBOccurCount/total + eps);%计算Y所蕴含的自信息（熵）
end

NMI = 2 * MI / (Hx+Hy);
end

