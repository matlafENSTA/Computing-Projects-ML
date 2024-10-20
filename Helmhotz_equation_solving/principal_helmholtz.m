%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRINCIPAL_HELMHOLTZ;                                                        %
%                                                                             %
% Une routine pour la mise en oeuvre des EF P1 Lagrange pour resoudre         %
% l'equation de Helmholtz stationnaire suivante                               %
%                                                                             %
% -\rho \omega^2  P - div(\mu \grad P) = S,   dans \Omega=\Omega_1 U \Omega_2 %
%                                                                             %
% ou S est la source, omega la pulsation,                                     %
% \rho constante, et \mu variable :                                           %
%                       \mu = | \mu_1,  dans \Omega_1                         %
%                             | \mu_2,  dans \Omega_2.                        %
%                                                                             %
% A l'equation de Helmholtz s'ajoute l'une des conditions suivantes :         %
%     1) une condition de Dirichlet :                                         %
%                         P = P_\Gamma,   sur le bord,                        %
%        ou P_\Gamma represente la pression exterieure ;                      %
%                                                                             %
%     2) une condition de Fourier :                                           %
%                        dP/dn + beta P = 0 sur le bord,                      %
%        ou beta est le coefficient de Fourier ;                              %
%                                                                             %
%     3) des conditions mixtes de Dirichlet-Fourier.                          %
%                                                                             %
% La routine calcule et affiche la solution, qui peut ensuite etre comparee   %
% a une solution de reference.                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------------- %
% Donnees du probleme %
% ------------------- %
nom_maillage = 'maillages/geom-2-couches.msh';
validation = 'non';           % Pour valider 'oui' ou 'non' les calculs

% Choix des conditions aux limites
CL = 'Dirichlet'; % 'Dirichlet', 'Neumann' ou 'Fourier'ou 'Mixtes'

% Donnees du probleme
rho = 1.0;        % Coefficient du probleme
omega = 5;      % Frequence
PP_Gamma = 1;   % Pression exterieure
alpha = -rho * omega^2;

if (strcmp(CL, 'Fourier') || strcmp(CL, 'fourier') ||...
    strcmp(CL, 'Mixtes')  || strcmp(CL, 'mixtes'))
    beta = -1i * omega;   % Condition de Fourier
    % beta = 0;           % Condition de Neumann
end

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
LL = zeros(Nbpt, 1);     % vecteur second membre

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

% Matrice de masse surfacique pour les conditions de Fourier %
% ---------------------------------------------------------- %
if (strcmp(CL, 'Fourier') || strcmp(CL, 'fourier') ||...
    strcmp(CL, 'Mixtes')  || strcmp(CL, 'mixtes'))
    % On trouve les aretes sur lesquelles la condition de Fourier est imposee
    idFourier = find(Refaretes == 2);   % Indices des aretes pour Fourier
    NbFourier = length(idFourier);      % Nombre d'aretes sur lesquelles Fourier est imposee
    SS = sparse(Nbpt, Nbpt);            % Matrice de masse surfacique

    % Assemblage sur les aretes de reference 2
    for l = 1:NbFourier
        % Matrice de masse surfacique elementaire
        % LA ROUTINE matS_elem.m DOIT ETRE COMPLETEE
        [Sel] = matS_elem(Coorneu(Numaretes(idFourier(l),1),:),...
                          Coorneu(Numaretes(idFourier(l),2),:), Refaretes(idFourier(l)));

        % On fait l'assemblage de la matrice de masse surfacique globale
        % A COMPLETER
    end % for l
end % for if

% ---------- %
% Matrice EF %
% ---------- %
AA = alpha*MM + KK;

if (strcmp(CL, 'Fourier') || strcmp(CL, 'fourier') ||...
    strcmp(CL, 'Mixtes')  || strcmp(CL, 'mixtes'))
  % On rajoute la contribution de la matrice de masse surfacique
  % dans le cas des condtions de 'Fourier' ou 'Mixtes'
  % COMPLETER
  % AA = ...
end

% ------------------------- %
% Calcul du second membre F %
% ------------------------- %
%  /!\ ATTENTION : f prend aussi omega en argument /!\
FF = f(Coorneu(:,1),Coorneu(:,2),omega);
LL = MM*FF;

% ------------------------------- %
% Pseudo-elimination et inversion %
% ------------------------------- %
if (strcmp(CL, 'Dirichlet') || strcmp(CL, 'dirichlet') ||...
    strcmp(CL, 'Mixtes')    || strcmp(CL, 'mixtes'))
  % Conditions de Dirichlet ou conditions mixtes
  [tilde_AA, tilde_LL] = elimine(AA, LL, Refneu);
  % tilde_AA ET tilde_LL SONT LA MATRICE EF ET LE VECTEUR SECOND MEMBRE
  % APRES PSEUDO_ELIMINATION
  UU = tilde_AA \ tilde_LL;
  PP = UU + PP_Gamma; % ??????????????
else
  % Conditions de Neumann ou de Fourier
  PP = AA \ LL;
end

if strcmp(validation, 'oui')
  % ------------------------------------------------------------------ %
  % Pour tracer la norme de l'inverse de A(omega) en fonction de omega %
  % ------------------------------------------------------------------ %
  %la fonction helmholtz_omega.m est dédiée à ce calcul.
end

% ---------- %
% Validation %
% ---------- %
if strcmp(validation, 'oui')
  % Solution de reference
  PP_exact = PP_Gamma + sin(3*pi*Coorneu(:,1)) .* sin(2*pi*Coorneu(:,2));
  % Calcul de l'erreur L2
  % A COMPLETER
  % Calcul de l erreur H1
  % A COMPLETER
  % attention de bien changer le terme source (dans FF)
end

% ------------- %
% Visualisation %
% ------------- %
if (strcmp(CL, 'Dirichlet') || strcmp(CL, 'dirichlet'))
  % La solution est reelle : on l'affiche donc dans une seule figure
  figure;
  affiche(real(PP), Numtri, Coorneu, sprintf('Dirichlet - %s', nom_maillage));
  set(gca, 'DataAspectRatio',[1 1 1]);
  %pour faire varier mu facilement dans l'affichage
  text(1.5, -0.3, sprintf('source radiale pour f, mu_1 = %d et mu_2 = %d',mu_1,mu_2), 'HorizontalAlignment', 'center');
end

if (strcmp(CL, 'Fourier') || strcmp(CL, 'fourier') ||...
    strcmp(CL, 'Mixtes')  || strcmp(CL, 'mixtes'))
  % La solution est COMPLEXE : on affiche donc ses parties reelle et
  % imaginaire dans des figures
  figure;
  affiche(real(PP), Numtri, Coorneu, sprintf('%s - Real(P) - %s', CL, nom_maillage));
  set(gca, 'DataAspectRatio',[1 1 1]);

  figure;
  affiche(imag(PP), Numtri, Coorneu, sprintf('%s - Imag(P) - %s', CL, nom_maillage));
  set(gca, 'DataAspectRatio',[1 1 1]);
end

% Affichage de la solution exacte et de la solution approchée pour comparer
if (strcmp(validation, 'oui'))
    figure;
    subplot(1,2,1);
    affiche(real(PP_exact), Numtri, Coorneu, sprintf('Interpolée de la solution exacte - %s - %s', CL, nom_maillage));
    set(gca, 'DataAspectRatio',[1 1 1]);
    subplot(1,2,2); %%%%%% attention, à enlever si changement de conditions aux limites !!!
    affiche(PP, Numtri, Coorneu, sprintf('Solution approchée - %s - %s', CL, nom_maillage));
    set(gca, 'DataAspectRatio',[1 1 1]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2023
