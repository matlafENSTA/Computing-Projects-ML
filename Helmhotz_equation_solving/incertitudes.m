function [logiL2,logiH1] = incertitudes(UU, MM, KK, Coorneu)

u = sin(3*pi*Coorneu(:,1)) .* sin(2*pi*Coorneu(:,2));%vecteur-valeurs de la
% solution exacte en tous les noeuds du maillage

%calcul de l'erreur L2
icL2 = (u-UU)'*MM*(u-UU);
logiL2 = log(icL2/(u'*MM*u))/2;
%calcul de l'erreur avec la semi-norme H1
icH1 = (u-UU)'*KK*(u-UU);
logiH1 = log(icH1/(u'*KK*u))/2;
