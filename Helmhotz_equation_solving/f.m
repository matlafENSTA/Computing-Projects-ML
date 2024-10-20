function val = f(x, y, omega)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F :
% Evaluation de la fonction second membre.
%
% SYNOPSIS val = f(x, y, omega)
%
% INPUT * x,y : les 2 coordonnees du point ou on veut evaluer la fonction.
%         omega : la pulsation
%
% OUTPUT - val: valeur de la fonction sur ce point.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choix de l'expression de la source
% A CHANGER EVENTUELLEMENT POUR LES SIMULATIONS
choix_source = 'radiale'; % 'sin-sin', 'constante', 'radiale', ou 'radiale'

if strcmp(choix_source,'validation')
  %val = ones(size(x));
  val = (13*pi^2-omega^2)*sin(3*pi*x).*sin(2*pi*y);
end

if strcmp(choix_source,'sin-sin')
  % Polynome trigonometrique
  val = sin(3*pi*x).*sin(2*pi*y);
end

if strcmp(choix_source,'constante')
  % Terme source constant
  val = 1;
end

if strcmp(choix_source,'radiale')
  % Source radiale et sinusoidale
  A = 5;           % Amplitude de la source
  x0 = 2.0;        % Abscisse du centre (2.0 pour Exos1-2, -0.5 pour Exo3)
  y0 = 0.5;        % Ordonnee du centre (0.5 pour Exos1-2,  0.0 pour Exo3)
  rr = sqrt((x - x0).^2 + (y - y0).^2);   % Variable radiale
  val = A * cos(omega*rr);                % Expression de la source
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2023
