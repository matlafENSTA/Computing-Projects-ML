function [Sel] = matS_elem(S1, S2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matS_elem :
% calcul la matrice de masse surfacique
%
% SYNOPSIS [Sel] = mat_elem_surface(S1, S2)
%
% INPUT * S1, S2 : les 2 coordonnees des 2 sommets de l'arete
%                      (vecteurs reels 1x2)
%
% OUTPUT - Sel matrice de masse surfacique elementaire pour le bord.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% preliminaires, pour faciliter la lecture:
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);

% On calcule la longueur de l'arete.
L = sqrt((x2 - x1)^2 + (y2 - y1)^2);

% Declaration et remplissage des matrice elementaires.
%---------------------------------------------------
Sel = zeros(2,2);
for i = 1:2
  for j = 1:2
    % A COMPLETER
    Sel(i, j) = ....
  end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2023
