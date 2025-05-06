%% Получение и визуализация эксперимента для клеток A549

load('eps_t.mat');
time = round(time - 0.999,12);
figure();
p = plot(time, strain,'k-o');
    xlabel('t', FontSize=14);
    ylabel('\epsilon(t)', FontSize=14);
    p(1).LineWidth = 0.7;
 
% применение различных методов сглаживания данных и визуализация
meanDef = {movmean(strain,3); movmean(strain,4); movmean(strain,5); movmean(strain,6)};
medianDef = {movmedian(strain,3); movmedian(strain,4); movmedian(strain,5); movmedian(strain,6)};

%  figure();
% p = plot(time, strain,'k-',time, meanDef{1}, 'b', time, meanDef{2}, 'r',time, meanDef{3}, 'g', time, meanDef{4}, 'm');
%     xlabel('t', FontSize=14);
%     ylabel('\epsilon(t)', FontSize=14);
%     p(1).LineWidth = 2;
%     p(2).LineWidth = 1.6;
%     p(3).LineWidth = 1.6;
%     p(4).LineWidth = 1.6;
%     p(5).LineWidth = 1.6;
%     grid on;
%     legend({'Original data','movmean with window = 3','movmean with window = 4','movmean with window = 5','movmean with window = 6'},'Location','southeast');

%  figure();
% p1 = plot(time, strain,'k-',time, medianDef{1}, 'b', time, medianDef{2}, 'r',time, medianDef{3}, 'g', time, medianDef{4}, 'm');
%     xlabel('t', FontSize=14);
%     ylabel('\epsilon(t)', FontSize=14);
%     p1(1).LineWidth = 2;
%     p1(2).LineWidth = 1.6;
%     p1(3).LineWidth = 1.6;
%     p1(4).LineWidth = 1.6;
%     p1(5).LineWidth = 1.6;
%     grid on;
%     legend({'Original data','movmedian with window = 3','movmedian with window = 4','movmedian with window = 5','movmedian with window = 6'},'Location','southeast');    

smoothDef = {smoothdata(strain,"sgolay"); smoothdata(strain,"lowess"); smoothdata(strain,"loess"); smoothdata(strain,"rlowess"); smoothdata(strain,"rloess")};
% фильтр Савитского-Голея, локальная регрессия с использованием взвешенных линейных наименьших квадратов и полиномиальной модели 1-й степени, 
% локальная регрессия ... полиномиальной модели 2-й степени,надёжная версия
% 'lowess' (меньший вес выбросам в регрессии), надёжная версия 'loess' (меньший вес выбросам в регрессии)
% smoothDef{1}(1) = 0; smoothDef{2}(1) = 0; smoothDef{3}(1) = 0; smoothDef{4}(1) = 0; smoothDef{5}(1) = 0;

NormEv =  [sqrt(sum((strain - smoothDef{1}).^2)); sqrt(sum((strain - smoothDef{2}).^2)); sqrt(sum((strain - smoothDef{3}).^2)); sqrt(sum((strain - smoothDef{4}).^2)); sqrt(sum((strain - smoothDef{5}).^2))];

 figure();
p2 = plot(time, strain,'k-', time, smoothDef{1}, 'b-', time, smoothDef{2}, 'r-',time, smoothDef{3}, 'g:', time, smoothDef{4}, 'm:', time, smoothDef{5}, ':');
    xlabel('t', FontSize=14);
    ylabel('\epsilon(t)', FontSize=14);
    p2(1).LineWidth = 2.5;
    p2(2).LineWidth = 2.5;
    p2(3).LineWidth = 2.5;
    p2(4).LineWidth = 2.5;
    p2(5).LineWidth = 2.5;
    p2(6).LineWidth = 2.5;
    grid on;
    legend({'Original data','Savitzky-Golay method', 'Local regression 1st DPM', 'Local regression 2nd DPM', 'Robust lowess', 'Robust loess'}, 'Location','southeast');

 figure();
 p3 = plot(time, strain,'k-', time, smoothDef{2}, 'r-');
    xlabel('t', FontSize=14);
    ylabel('\epsilon(t)', FontSize=14);
    p3(1).LineWidth = 2;
    p3(2).LineWidth = 2;
    grid on;
    legend({'Original data','Local regression 1st DPM'}, 'Location','southeast');

%% Определение частоты собственных колебаний

L = [50:1:280] * 1e-9;
E = [1, 2.5, 4] * 1e9;

f = FindOscillationFrequency(L, E);

figure();
p = plot(L*1e9, f(1,:)*1e-9, L*1e9, f(2,:)*1e-9, L*1e9, f(3,:)*1e-9); %зависимость длины филамента L от частоты колебаний при различных значениях E
    xlabel('L, nm', FontSize=12);
    ylabel('f, GHz',FontSize=12);
    p(1).LineWidth = 2;
    p(2).LineWidth = 2;
    p(3).LineWidth = 2;
    grid on;
legend({'E = 1 GPa','E = 2.5 GPa','E = 4 GPa'},'Location','northeast');

%% Построение кривых потенциала Psi

xi = 0.1; 
eRange = -0.3:0.01:0.8;
sigm0d = [-0.6 -0.1 0.3]';
F = Fpotential(sigm0d, eRange, xi);
figure();
pp = plot(eRange, F(1,:), eRange, F(2,:), eRange, F(end,:));
    pp(1).LineWidth = 2;
    pp(2).LineWidth = 2;
    pp(3).LineWidth = 2;
    xlim([-0.5 1]);
    legend({'sod = -0.6','sod = -0.1','sod = 0.3'},'Location','northwest');
dF = deltaF(sigm0d,xi);
nnu = nuCalculation(sigm0d,xi,Ge,tau0);

%% Идентификация параметров модели для напряжения sigm = 0.298 Па (200 мл/час)

close all; 
expDef = smoothDef{2};
t = time;
thet = 300 * 1.38 * 10^-23;
Gr = 23;
Ge = 0.4e9;
L = 5e-6; %длина сегмента
gam = 2.5 * 10e-22;
xi = 0.1; % xi должно быть = 0.1 (иначе кривые потенциала в диапазоне s0d [-0.6; -0.1; 0.3] не получаются качественно различными)
sigm = 0.298;
eps_init = expDef(1);

eps = FindDeformation(t, Gr, Ge, gam, L, thet, sigm, eps_init); %определение зависимости eps(t) с нач. значениями
figure();
p1 = plot(time, expDef,'ko-', t, eps,'r'); 
    xlabel('t, c', FontSize=12);
    ylabel(char(949),FontSize=14);
    p1(1).MarkerSize = 3;
    p1(1).LineWidth = 1;
    p1(1).MarkerFaceColor = 'k';
    p1(2).LineWidth = 2;
    grid on;
    legend({'Experimental data for \sigma = 0.298 Pa','Numerical results with identified parameters'},'Location','southeast');

% eps0 = eps - ((sigm - Gr * eps) / ( Ge));
% xi = 0.1;
% lamd = thet / (xi * gam);
% tau0 = 1/(FindOscillationFrequency(L, 3 * Ge));
% figure();
% p2 = plot(time, eps0);
%     p2(1).LineWidth = 2;
%     legend({'eps0'},'Location','southeast');
% 
% sigm0 = sigm - eps * Gr;
% figure();
% p3 = plot(time, sigm0);
%     p3(1).LineWidth = 2;
%     legend({'sigm0'},'Location','southeast');
% 
% sigm0d = sigm0 ./ (xi * lamd); % напряжения в безразмерном виде
% nu = nuCalculation(sigm0d, xi, Ge, tau0);
% figure();
% p4 = plot(time, nu);
%     p4(1).LineWidth = 2;
%     legend({'nu'},'Location','southeast');

x0 = [Gr Ge]; %решение оптимизационной задачи и идентификация параметров модели
options = optimset('MaxIter', 10^6,'MaxFunEvals',10^6);
[xSol, valOF] = fminsearch(@(x) FindObjectiveFunction(t, expDef, x(1), x(2), gam, L, thet, sigm, eps_init), x0, options);
numDef = FindDeformation(t, xSol(1), xSol(2), gam, L, thet, sigm, eps_init);

figure();
p1 = plot(time, expDef,'ko-', t, numDef,'r'); %сопоставление численных (при идентифицированных параметрах) и экспериментальных результатов
    xlabel('t, c', FontSize=12);
    ylabel(char(949),FontSize=14);
    p1(1).MarkerSize = 3;
    p1(1).LineWidth = 1;
    p1(1).MarkerFaceColor = 'k';
    p1(2).LineWidth = 2;
    xlim([0 42]);
    grid on;
    legend({'Experimental data for \sigma = 0.298 Pa','Numerical results with identified parameters'},'Location','southeast');

w_err = max((abs(expDef - numDef)./expDef) * 100);

%% Исследование устойчивости
% устойчивость по начальным данным

eps_true = numDef; %истинное решение (без возмущений)
Gr = xSol(1);
Ge = xSol(2);

n = 500; % число выборок для задания различных возмущений
nb = 5; % число возмущений
nb_step = 20; % шаг между возмущениями (в процентах)
% возмущения должны образовывать арифметическую прогрессию типа 1% 2% 3% 4% 5% => nb = 5, nb_step = 1 ; 5% 10% 15% =>nb = 3, nb_step = 5


bounds = eps_init * [ones(1,nb)-(1:nb) * nb_step * 0.01; ones(1,nb)+(1:nb) * nb_step * 0.01 ]'; %список границ возмущений н.у. 
ei_sampl = bounds(:,1) + (bounds(:,2) - bounds(:,1)) .* rand(1,n); % массив выборок н.у. (равномерное распределение на границах bounds)

figure(); %проверка выборочных данных (столбчатые диаграммы для различных возмущений н.у.)
tiledlayout(2,3);
for k = 1:5
    nexttile;
    histogram(ei_sampl(k,:),50);
    title(['Возмущение н.у. ', num2str(k), ' %']);
end

norm_list = zeros(5,n);
for k = 1:5
    for i = 1:n
        eps_pert = FindDeformation(t, Gr, Ge, gam, L, thet, sigm, ei_sampl(k, i)); %возмущенное решение
        norm_value = norm(eps_pert - eps_true) / norm(eps_true); %евклидова норма отклонения решений, записанная в относительных единицах (после деления)
        norm_list(k, i) = norm_value;
    end
end

figure(); %иллюстрация распределений норм отклонений решений при различных возмущениях 
tiledlayout(2,3);
for k = 1:5
    nexttile;
    histogram(norm_list(k,:), 50);
    title(['Отклонение при ', num2str(k), ' %']);
end

mean_list = [mean(norm_list(1,:)) mean(norm_list(2,:)) mean(norm_list(3,:)) mean(norm_list(4,:)) mean(norm_list(5,:))]; %мат. ожидание для каждого возмущения 1%...5%
std_list = [std(norm_list(1,:)) std(norm_list(2,:)) std(norm_list(3,:)) std(norm_list(4,:)) std(norm_list(5,:))]; %СКО для каждого возмущения 1%...5%

figure(); %построение графиков МО и СКО для различных возмущений 1%..5%
p1 = plot(1:5, mean_list, 'k-o');
    set(gca, 'XTick',1:5)
    xlabel('Δε₀, %', FontSize=12);
    ylabel('M',FontSize=14);
    p1.LineWidth = 2.5;
    p1.MarkerSize = 9;
    grid on;
figure();
p2 = plot(1:5, std_list, 'k-o');
    set(gca, 'XTick',1:5)
    xlabel('Δε₀, %', FontSize=12);
    ylabel('S',FontSize=14);
    p2.LineWidth = 2.5;
    p2.MarkerSize = 9;
    grid on;

%% локальные функции
function valOF = FindObjectiveFunction(tExp, eExp, Gr, Ge, gam, L, thet, sigm, epsInit) %найти значение целевой функции при известных экспериментальных деформациях yExp и векторе параметров xParam
eNum = FindDeformation(tExp, Gr, Ge, gam, L, thet, sigm, epsInit);
valOF = sqrt(sum(((eExp-eNum)./eExp).^2));
end

function def = FindDeformation(t, Gr, Ge, gam, L, thet, sigm, epsInit) %зависимость общей деформации(t) при постоянном напряжении: xP=[G0, G1, sigm, lamd, gam, thet, chi, tau0]
[~, eps] = ode15s(@(t, eps) odeF(eps, Gr, Ge, gam, L, thet, sigm), t, epsInit);
def = eps;
end

function dedt = odeF(eps, Gr, Ge, gam, L, thet, sigm) %функция производной для решения ДУ в FindDeformation
xi = 0.1;
lamd = thet / (xi * gam);
tau0 = 1/(FindOscillationFrequency(L, 3 * Ge));
eps0 = eps - ((sigm - Gr * eps) / ( Ge));
sigm0 = sigm - eps * Gr;
ksiVal = ksiDependence(eps0);
sigm0d = sigm0 ./ (xi * lamd); % напряжения в безразмерном виде
nu = nuCalculation(sigm0d, xi, Ge, tau0);
dedt = (Ge / (nu * (Ge + Gr))) * (sigm - Gr * eps - (thet/(gam)) * ksiVal + lamd * eps0);
% g1 = sigm - Gr * eps;
% g2 = - (thet/(gam)) * ksiVal + lamd * eps0;
end

function f = FindOscillationFrequency(L, E)
% L - длина актинового филамента, E - модуль Юнга
beta = 120.9;
r = 4e-9; % радиус поперечного сечения
roA = 2.53e-14;
I = pi * r^4 / 4;
f = (beta./(2 * pi * L.^2)) .* sqrt(E' .* I / (roA));
end

