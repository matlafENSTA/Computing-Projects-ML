function val = mu_1(x,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MU_1 :
% Evaluation de la fonction mu_1.
%
% SYNOPSIS val = mu_1(x,y)
%
% INPUT * x,y : les 2 coordonnees du point ou on veut evaluer la fonction.
%
% OUTPUT - val: valeur de la fonction sur ce point.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choix de l'expression du coefficient mu_1
% A CHANGER EVENTUELLEMENT POUR LES SIMULATIONS
choix_mu1 = 'constante';   % 'constante' ou 'variable'


if strcmp(choix_mu1, 'constante')
  % Coefficient constant
  val = 1;
end

if strcmp(choix_mu1,'variable')
  % Interface penchee / mu1 variable et stratifie le long de l'interface
  % Source dans Omega2
  p = 3;
  mu2 = 1;    % Coefficient dans Omega2
  C = 3*mu2;  % Valeur sur l'interface; si C = mu2, alors le saut est nul
  A = 4*C/5;  % Amplitude
  alpha = (2*p+1)*pi;
  a = 1/4.;   % Pente de l'interface
  val = C + A* cos(alpha*(y-a*x));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2023
