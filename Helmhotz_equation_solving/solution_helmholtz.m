function [UU, MM, KK] = solution_helmholtz(nom_maillage, omega, beta, visualisation)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % SOLUTION_HELMHOLTZ;
  %
  % Calcul d'une solution approchee par elements finis P1 de l'equation de
  % Helmholtz stationnaire, avec conditions mixtes de Dirichlet homogene sur
  % une partie du bord, et de Fourier homogene sur l'autre partie.
  %
  % | -\omega^2  U - div(\mu \grad U) = f,   dans \Omega
  % |                              U  = 0,   sur  \Gamma_D
  % |                dU/dn + \beta U  = 0,   sur  \Gamma_F,
  %
  % ou f est la source, omega la pulsation, \mu variable, \beta le coefficient
  % de transmission, et ou (\Gamma_D) U (\Gamma_G) = \Gamma, le bord de \Omega.
  %
  % SYNOPSIS [UU, MM, KK] = solution_helmholtz(nom_maillage, omega, beta, visualisation)
  %
  % INPUT  (1) nom_maillage, le nom du maillage sur lequel est resolue l'EDP,
  %        (2) omega, la pulsation,
  %        (3) beta, le coefficient de transmission,
  %        (4) visualisation, booleen qui indique si on veut voir ou non la
  %            solution et ses derivees.
  %
  % OUTPUT (1) UU, vecteur contenant les valeurs de la solution approchee en
  %            les points du maillage (de taille Nbpt x 1),
  %        (2) MM, matrice de masse,
  %        (3) KK, matrice de rigidite.
  %
  % NOTE   (1) La condition de Dirichlet est imposee sur les aretes de
  %            reference 1, et celle de Fourier, sur les aretes de reference 2.
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % --------------------- %
  % Parametres par defaut %
  % --------------------- %
  if (nargin < 4)
    visualisation = 0; % Ne pas afficher de solution par defaut
  end
  if (nargin < 3)
    beta = 0; % Conditions mixtes de Dirichlet-Neumann par defaut
  end
  if (nargin < 2)
    omega = 1.0;
  end

  % -------------------------------- %
  % Lecture du maillage et affichage %
  % -------------------------------- %
  [Nbpt, Nbtri, Coorneu, Refneu, Numtri, ...
   Reftri, Nbaretes, Numaretes, Refaretes] = lecture_msh(nom_maillage); %#ok

  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Debut de la partie a completer %
  % ---------------------- %
  % Calcul des matrices EF %
  % ---------------------- %
  % Declarations %
  MM = sparse(Nbpt, Nbpt); % matrice de masse
  KK = sparse(Nbpt, Nbpt); % matrice de rigidite

  % Boucle sur les triangles %
  % ------------------------ %
  for l = 1:Nbtri
      % Matrice de masse elementaire
      [Mel] = matM_elem(Coorneu(Numtri(l,1),:),Coorneu(Numtri(l,2),:), ...
                        Coorneu(Numtri(l,3),:));

      % Matrice de rigidite elementaire
      % LA ROUTINE matK_mu_elem.m DOIT ETRE MODIFIEE
      [Kel] = matK_mu_elem(Coorneu(Numtri(l,1),:), Coorneu(Numtri(l,2),:), ...
    		                   Coorneu(Numtri(l,3),:), Reftri(l));

      % On fait l'assemblage de la matrice globale
      % A COMPLETER
  end % for l

  % Matrice de masse surfacique pour les conditions de Fourier %
  % ---------------------------------------------------------- %
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

  % ---------- %
  % Matrice EF %
  % ---------- %
  alpha = -omega * omega;
  % On rajoute la contribution de la matrice de masse surfacique
  % dans le cas des condtions de 'Fourier' ou 'Mixtes'
  % COMPLETER
  AA = ...

  % ------------------------- %
  % Calcul du second membre F %
  % ------------------------- %
  % A COMPLETER EN UTILISANT LA ROUTINE f.m
  %  /!\ ATTENTION : f prend aussi omega en argument /!\
  FF = f(...);
  LL = ...;

  % ------------------------------------------------------------- %
  % Pseudo-elimination sur les aretes de reference 1 et inversion %
  % ------------------------------------------------------------- %
  % tilde_AA ET tilde_LL SONT LA MATRICE EF ET LE VECTEUR SECOND MEMBRE
  % APRES PSEUDO_ELIMINATION
  % ECRIRE LA ROUTINE elimine.m ET INSERER L'APPEL A CETTE ROUTINE
  % A UN ENDROIT APPROPRIE
  UU = tilde_AA \ tilde_LL;
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fin de la partie a completer %

  if (visualisation)
    % Affichage de la solution %
    % ------------------------ %
    figure;
    subplot(1, 2, 1);
    affiche(real(UU), Numtri, Coorneu);
    xlim([-1, 1]); ylim([-1, 1]);
    title('Partie reelle de la solution');
    set(gca, 'DataAspectRatio',[1 1 1]);

    subplot(1, 2, 2);
    affiche(imag(UU), Numtri, Coorneu);
    xlim([-1, 1]); ylim([-1, 1]);
    title('Partie imaginaire de la solution');
    set(gca, 'DataAspectRatio',[1 1 1]);

    % Approximation et affichage des derivees partielles %
    % -------------------------------------------------- %
    DD1 = sparse(Nbpt, Nbpt); %
    DD2 = sparse(Nbpt, Nbpt); %
    for l = 1:Nbtri
        % Matrice EF elementaire
        [Del1, Del2] = matD_elem(Coorneu(Numtri(l,1),:),Coorneu(Numtri(l,2),:), ...
                                 Coorneu(Numtri(l,3),:));

        % Assemblage
        DD1(Numtri(l,:), Numtri(l,:)) = DD1(Numtri(l,:), Numtri(l,:)) + Del1; %#ok
        DD2(Numtri(l,:), Numtri(l,:)) = DD2(Numtri(l,:), Numtri(l,:)) + Del2; %#ok
    end % for l

    % Affichage des parties reelles
    figure;
    subplot(1, 2, 1);
    affiche(real(MM\(DD1*UU)), Numtri, Coorneu);
    xlim([-1, 1]); ylim([-1, 1]);
    title('Partie reelle de la derivee par rapport a x');
    set(gca, 'DataAspectRatio',[1 1 1]);

    subplot(1, 2, 2);
    affiche(real(MM\(DD2*UU)), Numtri, Coorneu);
    xlim([-1, 1]); ylim([-1, 1]);
    title('Partie reelle de la derivee par rapport a y');
    set(gca, 'DataAspectRatio',[1 1 1]);
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2023
