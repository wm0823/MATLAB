function [alphas,offset] = my_smo(data, targetLabels, boxConstraints, smoOptions)

% 第一个alpha 的选择方法为 搜索第一个不满足KKT条件的
% 第二个alpha 选择分两步，首先选择|e1 - e2|，若满足条件则保留，否则随机选择
%最大迭代次数为smoOptions.MaxIter，或者当检查一遍没有变量进一步优化，则结束

%初始化
tolKKT = smoOptions.TolKKT;
nPoints = length(targetLabels);
alphas = zeros(nPoints, 1);
offset = 0;
gx = zeros(1,nPoints);
ei = zeros(1,nPoints);
itCount =0;
maxIter = smoOptions.MaxIter;
% 迭代
while (itCount < maxIter)
    flag = 0;
    for outai = 1: nPoints
      
        for i = 1 : nPoints
            gx(i) = (alphas .*  targetLabels)' * data * data(i, :)' + offset;
            ei(i) = gx(i) - targetLabels(i);
        end
        
        %选择第一个alpha
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
             %选择第二个alpha 
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
            
           
            %更新alpha2
           nalpha2 = alphas(indexa2) + targetLabels(indexa2)*(ei1-ei2)/eta;
            if nalpha2 < Low
                nalpha2 = Low;
            elseif nalpha2 > High
                nalpha2 = High;
            end
            
            if abs(nalpha2 - alphas(indexa2)) < 1e-3
                continue;
            end
            %更新alpha1的值
            nalpha1 = alphas(indexa1) + targetLabels(indexa1)*targetLabels(indexa2)*(alphas(indexa2)-nalpha2);
            %更新offset
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
disp('迭代次数');
disp(itCount);
end
%检查KKT条件
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
%随机选择
function index = selectJrand(i,m)
    index = i;
    while(index == i)
        index = randi([1,m],1,1);  
    end
end