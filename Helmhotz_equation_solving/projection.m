function UU2 = projection(UU1, nom_maillage_1, nom_maillage_2)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PROJECTION; Projection d'un vecteur sur un autre maillage.
  %
  % SYNOPSIS UU2 = projection(UU1, nom_maillage_1, nom_maillage_2)
  %
  % INPUT  (1) UU1, vecteur sur le maillage n°1 (de taille Nbpt1 x 1),
  %        (2) nom_maillage_1, maillage n°1,
  %        (3) nom_maillage_2, maillage n°2.
  %
  % OUTPUT (1) UU2, vecteur correspondant a la solution du vecteur UU1 projetee
  %            sur le maillage n°2 (de taille Nbpt2 x 1).
  %
  % NOTE   (1) nom_maillage_1 et nom_maillage_2 doivent correspondre au meme
  %            domaine.
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  [~, ~, Coorneu1, ~, Numtri1] = lecture_msh(nom_maillage_1);
  [Nbpt2, ~, Coorneu2, ~, ~]   = lecture_msh(nom_maillage_2);
  UU2 = zeros(Nbpt2, 1);
  tolerance = 1e-05;

  for i = 1:Nbpt2
    % On veut calculer UU2(i), qui correspond a la fonction de UU1 evaluee en
    % le i-eme noeud du maillage n°2.
    S2 = Coorneu2(i, :)';

    % Si le i-eme noeud du maillage n°2 est aussi un noeud du maillage n°1, il
    % suffit de trouver son indice dans ce maillage
    indiceS2maillage1 = find(abs(Coorneu1(:, 1) - S2(1)) + ...
                             abs(Coorneu1(:, 2) - S2(2)) < tolerance);

    if (length(indiceS2maillage1) == 1)
      % +-----+---------------------------------------- %
      % | /!\ |  Dans le TP, les maillages sont generes %
      % | /!\ |  de telle sorte qu'on soit dans ce cas  %
      % +-----+---------------------------------------- %
      UU2(i) = UU1(indiceS2maillage1);
      % ----------------------------------------------- %
    else
      fprintf(['Le point (%f, %f) n''est pas l''un des sommets du maillage 1',...
               'On utilise alors de l''interpolation\n'], S2(1), S2(2));
      warning(['Dans le TP ANN201, vous n''etes pas censes arriver a cet endroit.\n',...
               'Parlez-en a votre charge de TP !']);

      B = [P; 1];
      alpha = -inf;
      idtri = 1;

      while (true)
        % Noeuds du triangle
        S1tri = Coorneu1(Numtri1(idtri, 1), :)';
        S2tri = Coorneu1(Numtri1(idtri, 2), :)';
        S3tri = Coorneu1(Numtri1(idtri, 3), :)';

        % On dilate un peu le triangle courant pour rajouter une petite tolerance
        P1 = (1 + 2*tolerance)*S1tri - tolerance*(S2tri + S3tri);
        P2 = (1 + 2*tolerance)*S2tri - tolerance*(S3tri + S1tri);
        P3 = (1 + 2*tolerance)*S3tri - tolerance*(S1tri + S2tri);

        A = [P1, P2, P3; ones(1, 3)];
        alpha = A \ B;    % Coordonnees barycentriques du point dans chaque triangle

        if (alpha >= -tolerance && alpha <= 1 + tolerance)
          break;
        elseif idtri >= size(Numtri1, 1)
          error(['Le point (%f, %f) ne semble pas appartenir au domaine. ',...
                 'Augmenter la tolerance peut-etre'], S2(1), S1(2));
        else
          idtri = idtri + 1;
        end
      end

      UU2(i) = alpha.' * UU1(Numtri1(idtri, :));
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2023
