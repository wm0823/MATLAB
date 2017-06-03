function [alphas,offset] = my_smo(data, targetLabels, boxConstraints, smoOptions)

% ��һ��alpha ��ѡ�񷽷�Ϊ ������һ��������KKT������
% �ڶ���alpha ѡ�������������ѡ��|e1 - e2|�������������������������ѡ��
%����������ΪsmoOptions.MaxIter�����ߵ����һ��û�б�����һ���Ż��������

%��ʼ��
tolKKT = smoOptions.TolKKT;
nPoints = length(targetLabels);
alphas = zeros(nPoints, 1);
offset = 0;
gx = zeros(1,nPoints);
ei = zeros(1,nPoints);
itCount =0;
maxIter = smoOptions.MaxIter;
% ����
while (itCount < maxIter)
    flag = 0;
    for outai = 1: nPoints
      
        for i = 1 : nPoints
            gx(i) = (alphas .*  targetLabels)' * data * data(i, :)' + offset;
            ei(i) = gx(i) - targetLabels(i);
        end
        
        %ѡ���һ��alpha
        while outai < nPoints
            if checkKKT(alphas(outai), targetLabels(outai), boxConstraints(outai),gx(outai), tolKKT) == 0
                break;
            else
                outai = outai + 1;
            end
        end
        gx1 = gx(outai);
        ei1 = ei(outai);
     
        if checkKKT(alphas(outai), targetLabels(outai), boxConstraints(outai),gx1, tolKKT) == 0
             indexa1 = outai;
             %ѡ��ڶ���alpha 
             if ei(indexa1) > 0
                 [~, indexa2] = min(ei);
             else
                 [~, indexa2] = max(ei);
             end
            if indexa2 ~= indexa1 && alphas(indexa2) ~= 0 && alphas(indexa2) ~= boxConstraints(indexa2)
            else
                indexa2 = selectJrand(outai,nPoints);
                if indexa2 == indexa1
                    indexa2 = selectJrand(outai,nPoints);
                end
            end
            gx2 = gx(indexa2);
            ei2 = ei(indexa2);
            if targetLabels(indexa1)*targetLabels(indexa2) == 1
                Low = max(0, alphas(indexa2) + alphas(indexa1) - boxConstraints(indexa1));
                High = min(boxConstraints(indexa1), alphas(indexa2) + alphas(indexa1));
            else
                Low = max(0, alphas(indexa2) - alphas(indexa1));
                High = min(boxConstraints(indexa1), boxConstraints(indexa1) + alphas(indexa2) - alphas(indexa1));
            end
            
            if Low == High
                continue;
            end
            eta = data(indexa1, :)*data(indexa1,:)' + data(indexa2, :)*data(indexa2,:)' - 2*data(indexa1, :)*data(indexa2,:)';
           
            if eta <= 0
                continue;
            end
            
           
            %����alpha2
           nalpha2 = alphas(indexa2) + targetLabels(indexa2)*(ei1-ei2)/eta;
            if nalpha2 < Low
                nalpha2 = Low;
            elseif nalpha2 > High
                nalpha2 = High;
            end
            
            if abs(nalpha2 - alphas(indexa2)) < 1e-3
                continue;
            end
            %����alpha1��ֵ
            nalpha1 = alphas(indexa1) + targetLabels(indexa1)*targetLabels(indexa2)*(alphas(indexa2)-nalpha2);
            %����offset
            b1 = -ei1-targetLabels(indexa1)*(data(indexa1,:)*data(indexa1,:)')*(nalpha1-alphas(indexa1))- targetLabels(indexa2)*(data(indexa2,:)*data(indexa1,:)')*(nalpha2-alphas(indexa2))+offset;
            b2 = -ei2-targetLabels(indexa1)*(data(indexa1,:)*data(indexa2,:)')*(nalpha1-alphas(indexa1))- targetLabels(indexa2)*(data(indexa2,:)*data(indexa2,:)')*(nalpha2-alphas(indexa2))+offset;
            alphas(indexa2) = nalpha2;
            alphas(indexa1) = nalpha1;
            if ((alphas(indexa1)>0) && (alphas(indexa1)<boxConstraints(indexa1)))
                offset = b1;
            elseif ((alphas(indexa2)>0) && (alphas(indexa2)<boxConstraints(indexa2)))
                offset = b2;
            else
                offset = (b1+b2)/2.0;
            end
            flag = flag + 1;
        end
        
    end
    if flag == 0
        break;
    end

itCount = itCount +1;
end
disp('��������');
disp(itCount);
end
%���KKT����
function flag = checkKKT(Alpha, targetL, boxC, gxi, tol)
    flag = 0;
    if (Alpha == 0) && (targetL * gxi  >= 1 - tol)
        flag = 1;
    elseif (0 < Alpha) && (Alpha < boxC) && (abs(targetL * gxi - 1) == tol)
        flag = 1;
    elseif (Alpha == boxC) && (targetL * gxi  <= 1 + tol)
        flag = 1;
    end
end
%���ѡ��
function index = selectJrand(i,m)
    index = i;
    while(index == i)
        index = randi([1,m],1,1);  
    end
end