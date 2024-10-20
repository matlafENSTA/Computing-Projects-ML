function [Kel] = matK_mu_elem(S1, S2, S3, Reftri)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matK_mu_elem :
% Calcule la matrice de raideur elementaire en P1 lagrange avec un
% coefficient variable.
%
% SYNOPSIS [Kel] = matK_mu_elem(S1, S2, S3, Reftri)
%
% INPUT * S1, S2, S3 : les 2 coordonnees des 3 sommets du triangle
%                      (vecteurs reels 1x2)
%         Reftri : reference du triangle (1 si appartient a Omega_1
%                     et 2 si dans Omega_2)
%
% OUTPUT - Kel matrice de raideur elementaire (matrice 3x3)
%
% NOTE (1) le calcul est exact (pas de condensation de masse)
%      (2) calcul direct a partir des formules donnees par
%          les coordonnees barycentriques
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% preliminaires, pour faciliter la lecture:
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

% les 3 normales a l'arete opposees (de la longueur de l'arete)
norm = zeros(3, 2);
norm(1, :) = [y2-y3, x3-x2];
norm(2, :) = [y3-y1, x1-x3];
norm(3, :) = [y1-y2, x2-x1];

% D est, au signe pres, deux fois l'aire du triangle
D = abs((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1)); % = |det(B_l)|
if (abs(D) <= eps)
  error('l aire d un triangle est nulle!!!');
end
% calcul de B_l et S_l, matrices définissant F_l qui passe du triangle de
% référence au triangle Reftri.
B_l = [x2-x1,x3-x1;y2-y1,y3-y1];
S_l = [x1;y1];                        %colonne

%pour la formule de quadrature : 
w0 = 1/6;
s0 = 1/6;
s1 = 2/3;
%calcul des coordonnées des points images de (s0,s1),(s0,s0),(s1,s0) par
%F_l (fonction qui transforme le triangle de référence en T_l)
t00 = B_l*[s0;s0]+S_l;
t01 = B_l*[s0;s1]+S_l;
t10 = B_l*[s1;s0]+S_l;

%gradients des fonctions wi définies sur le triangle de référence:
grad_w1 = [-1;-1];
grad_w2 = [1;0];
grad_w3 = [0;1];
gw = [grad_w1,grad_w2,grad_w3];

% calcul de la matrice de raideur
% -------------------------------
% Attention a prendre en compte le coefficient variable mu et le sous-domaine
% Omega
Kel = zeros(3,3);
for i=1:3
  for j=1:3
    c = w0*abs(det(B_l))*(B_l'\gw(:,i))'*(B_l'\gw(:,j)); 
    if Reftri == 1
        Kel(i,j) = c*(mu_1(t00(1),t00(2)) + mu_1(t01(1),t01(2)) + mu_1(t10(1),t10(2)));
    elseif Reftri == 2
        Kel(i,j) = c*(mu_2(t00(1),t00(2)) + mu_2(t01(1),t01(2)) + mu_2(t10(1),t10(2)));
    end % j
  end % i
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2023
