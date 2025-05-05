function dF = deltaF(sod, xi) %выводит величину барьера потенциала F

epsoRange = -0.3:0.01:0.8;
F = Fpotential(sod,epsoRange,xi); %при получении массива sigm0 F предсталяет собой матрицу 
% и необходимо обрабатывать каждую строку (массив eRange с определенным значением sigm0)
dF = zeros(1,length(sod));
for i=1:length(sod)
    F1 = F(i,:);
    A1 = islocalmax(F1);
    A2 = islocalmin(F1);
    
    if(numel(epsoRange(A1)) + numel(epsoRange(A2)) == 0)
        if(F1(1) > F1(end))
            dF(i) = 0;
        else
            dF(i) = F1(end) - F1(1);
        end
    else
        if(numel(epsoRange(A1)) + numel(epsoRange(A2)) > 1) %если есть хотя бы 1 лок. минимум
            MaxInd = find(A1, 1); %выбираем самый левый максимум
            MinInd = find(A2, 1); %выбираем самый левый минимум
            if(epsoRange(MaxInd) > epsoRange(MinInd))
                dF(i) = F1(MaxInd) - F1(MinInd);
            else 
                disp("Error: Xmax < Xmin");
            end
        else
            if(F1(1) > F1(end))
                dF(i) = 0;
            else
                dF(i) = F1(find(epsoRange == 0.4)) - F1(A2);
            end
        end
    end
end
end


