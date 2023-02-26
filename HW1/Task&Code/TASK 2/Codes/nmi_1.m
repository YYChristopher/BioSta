function NMI = nmi_1( A, B ) 

if length( A ) ~= length( B)
    error('length( A ) must == length( B)');%�ж�ǰ�����������������ݵ���Ӧ��ͬ������ͬ��ֱ�ӽ���
end
total = length(A);
A_ids = unique(A);
B_ids = unique(B);%ȡ��A��B��һ����Ԫ�أ����ظ���

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
         
         MI = MI + pxy*log2(pxy/(px*py)+eps); % eps : ����ֵ����С�������˴��ȼ���MIֵ

    end
end

% Normalized Mutual information
Hx = 0; % Entropies
for idA = A_ids
    idAOccurCount = length( find( A == idA ) );
    Hx = Hx - (idAOccurCount/total) * log2(idAOccurCount/total + eps);%����X���̺�������Ϣ���أ�
end
Hy = 0; % Entropies
for idB = B_ids
    idBOccurCount = length( find( B == idB ) );
    Hy = Hy - (idBOccurCount/total) * log2(idBOccurCount/total + eps);%����Y���̺�������Ϣ���أ�
end

NMI = 2 * MI / (Hx+Hy);
end

