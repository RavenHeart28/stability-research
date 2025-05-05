function Fp = Fpotential(sod, eRange, xi) %выводит качественную зависимость, смещенную вблизи нуля
    Fp = real(-1.7*log(0.99-eRange) - 114161*log(8822.34-eRange) - 0.9*log(0.5+eRange) - 36.34*log(2.84+eRange)) - (eRange.^2)/(2*xi) - sod.* eRange;
    means = mean(Fp, 2); % найти среднее в каждой строке (на выходе столбец средних значений строк)
    Fp = Fp - means;
end

