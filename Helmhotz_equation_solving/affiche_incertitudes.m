function [] = affiche_incertitudes(omega)
%les valeurs de hi, logiL2i et logiH1i (logarithmes des erreurs pour chaque
% valeur de h et omega) ont été changées manuellement et classées dans 
% l'ordre décroissant de h. 

%omega = 1
%h1 = 0.8;
% logiL21 = -0.0420;
% logiH11 = -0.1209;
% h2 = 0.3;
% logiL22 = -0.9590;
% logiH12 = -1.0082;
% h3 = 0.1;
% logiL23 = -2.5690;
% logiH13 = -2.5235;
% h4 = 0.03;
% logiL24 = -4.9797;
% logiH14 = -4.8753;

%omega = 4
% h1 = 0.8;
% logiL21 = 0.5306;
% logiH11 = 0.2816;
% h2 = 0.3;
% logiL22 = -0.8012;
% logiH12 = -0.9095;
% h3 = 0.08;
% logiL23 = -2.9142;
% logiH13 = -2.9020;
% h4 = 0.03;
% logiL24 = -4.8537;
% logiH14 = -4.7719;

%omega = 4.44
% h1 = 0.9;
% logiL21 = 1.0331;
% logiH11 = 0.6358;
% h2 = 0.3;
% logiL22 = -0.4993;
% logiH12 = -0.8362;
% h3 = 0.08;
% logiL23 = -1.1224;
% logiH13 = -1.9978;
% h4 = 0.02;
% logiL24 = -5.0475;
% logiH14 = -5.4054;

%omega = 12.8 : sur un pic
h1 = 0.9;
logiL21 = -0.4786;
logiH11 = -0.4859;
h2 = 0.2;
logiL22 = 2.0818;
logiH12 = 2.0798;
h3 = 0.07;
logiL23 = -1.9182;
logiH13 = -1.9121;
h4 = 0.01;
logiL24 = -5.8467;
logiH14 = -5.8342;

%classement des valeurs
L2 = [logiL21,logiL22,logiL23,logiL24];
H1 = [logiH11,logiH12,logiH13,logiH14];
lnh = [log(1/h1),log(1/h2),log(1/h3),log(1/h4)];

%fonctions d'affichage
scatter(lnh,L2, 'r', 'filled', 'DisplayName', 'log(||u-uh||L2/||u||L2)');
hold on;
scatter(lnh,H1, 'b', 'filled', 'DisplayName', 'log(|u-uh|H1/|u|H1)');
xlabel('log(1/h)');
ylabel('différentes incertitudes');
hold on;

% Insertion d'une courbe de tendance pour l'erreur L2
coeffsL2 = polyfit(lnh, L2, 1);
% Création d'une fonction polynomiale à partir des coefficients
p = polyval(coeffsL2, lnh);
% Tracé des données expérimentales et de la courbe de tendance
plot(lnh, L2, 'o', lnh, p, '-');
title(sprintf('erreurs L2 et H1 pour omega = %f',omega));
grid on;