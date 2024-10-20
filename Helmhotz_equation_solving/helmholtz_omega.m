%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELMHOLTZ OMEGA VARIABLE (q1.16);                                           %
% ------------------- %
% Donnees du probleme %
% ------------------- %
nom_maillage = 'maillages/geomCarre.msh';
validation = 'oui';           % Pour valider 'oui' ou 'non' les calculs

% Choix des conditions aux limites
CL = 'Dirichlet'; % 'Dirichlet', 'Neumann' ou 'Fourier'ou 'Mixtes'

if strcmp(validation,'oui')
    mu_1 = 1;
    mu_2 = 1;
    PP_Gamma = 0;
end

% -------------------------------- %
% Lecture du maillage et affichage %
% -------------------------------- %
[Nbpt, Nbtri, Coorneu, Refneu, Numtri, ...
 Reftri, Nbaretes, Numaretes, Refaretes] = lecture_msh(nom_maillage);

% ---------------------- %
% Calcul des matrices EF %
% ---------------------- %
% Declarations %
% ------------ %
MM = sparse(Nbpt, Nbpt); % matrice de masse
KK = sparse(Nbpt, Nbpt); % matrice de rigidite

% Donnees du probleme
rho = 1.0;        % Coefficient du probleme
PP_Gamma = 0.0;   % Pression exterieure
OMEGA = [0.15:0.15:15];
norm_AA = zeros(size(OMEGA));
q = 1; %compteur pour l'indice de norm_AA
for omega = OMEGA % Frequence
    alpha = -rho * omega^2;
    
    % Boucle sur les triangles %
    % ------------------------ %
    for l = 1:Nbtri
        % Matrice de masse elementaire
        Mel = matM_elem(Coorneu(Numtri(l,1),:),Coorneu(Numtri(l,2),:),Coorneu(Numtri(l,3),:));
        % Matrice de rigidite elementaire
        Kel = matK_mu_elem(Coorneu(Numtri(l,1),:), Coorneu(Numtri(l,2),:),Coorneu(Numtri(l,3),:), Reftri(l));
      for i=1:3 % On fait l'assemblage de la matrice globale
          I = Numtri(l,i);
          for j=1:3
              J=Numtri(l,j);
              MM(I,J) = MM(I,J) + Mel(i,j);
              KK(I,J) = KK(I,J) + Kel(i,j);
          end %for j
      end %for i
    end % for l
    
    % ---------- %
    % Matrice EF %
    % ---------- %
    AA = alpha*MM + KK;
    
    % ------------------------------- %
    % Pseudo-elimination et inversion %
    % ------------------------------- %
    if (strcmp(CL, 'Dirichlet') || strcmp(CL, 'dirichlet') ||...
        strcmp(CL, 'Mixtes')    || strcmp(CL, 'mixtes'))
      % Conditions de Dirichlet ou conditions mixtes
      [tilde_AA, tilde_LL] = elimine(AA, LL, Refneu);
      % tilde_AA ET tilde_LL SONT LA MATRICE EF ET LE VECTEUR SECOND MEMBRE
      % APRES PSEUDO_ELIMINATION
    
    if strcmp(validation, 'oui')
      % ------------------------------------------------------------------ %
      % Pour tracer la norme de l'inverse de A(omega) en fonction de omega %
      % ------------------------------------------------------------------ %
      norm_AA(q)= sqrt(trace(inv(tilde_AA)'*inv(tilde_AA)));
      q = q + 1;
    end
    end
end
% ------------- %
% Visualisation %
% ------------- %
plot(OMEGA,norm_AA);
xlabel('omega');
ylabel('norm(inv(A))');
title('Comportement par rapport à omega pour un maillage avec h = 0,1');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2023
