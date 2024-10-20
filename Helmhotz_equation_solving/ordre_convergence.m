%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ORDRE_CONVERGENCE;
%
% Routine pour calculer les erreurs L2 et semi-norme H1 pour differents
% maillages, et pour les tracer en fonction du pas de maillage.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---------- %
% Parametres %
% ---------- %
omega = 1.0;                      % Frequence
beta = 1i*omega;                  % Coefficient de transmission
prefixe_maillage = 'maillages/geomPacman-030-';    % prefixe du maillage
Npts = 4;                         % Nombre de points de simulation
h = 0.5.^(0:Npts-1)';             % Pas de maillage
ErreursL2rel = zeros(Npts, 1);    % Erreur en norme L2
ErreursH1rel = zeros(Npts, 1);    % Erreur en semi-norme H1

% -------------------------------------------------- %
% Calcul de la solution de reference                 %
% La solution de reference est obtenue               %
% en utilisant la methode sur un maillage tres fin.  %
% -------------------------------------------------- %
fprintf('Calcul de la solution de reference\n');
nom_maillage_ref = [prefixe_maillage, '0.msh'];
UU_reference = solution_helmholtz(nom_maillage_ref, omega, beta, 1);

% ---------------------------------------------- %
% Calcul des solutions approchees et comparaison %
% ---------------------------------------------- %
for i = 1:Npts
  % Calcul de la solution
  fprintf('Calcul de la solution pour h = %f\n', h(i));
  nom_maillage = [prefixe_maillage, num2str(i), '.msh'];
  [UU, MM, KK] = solution_helmholtz(nom_maillage, omega, beta, 0);

  % Pour comparer solution approchee et solution de reference, il faut
  % projeter la solution de reference sur le maillage de la solution approchee
  UUref = projection(UU_reference, nom_maillage_ref, nom_maillage);

  % Calcul de l'erreur entre solution approchee et solution de reference
  ErreursL2rel(i) = sqrt(((UU - UUref)' * MM * (UU - UUref)) / (UUref' * MM * UUref));
  ErreursH1rel(i) = sqrt(((UU - UUref)' * KK * (UU - UUref)) / (UUref' * KK * UUref));
end

% --------------------------------------------------------- %
% Affichage des erreurs et calcul des ordres de convergence %
% --------------------------------------------------------- %
ErreursL2rel = real(ErreursL2rel);
ErreursH1rel = real(ErreursH1rel);
figure;
loglog(h, ErreursL2rel, 'ro'); hold on;
loglog(h, ErreursH1rel, 'bx');
xlabel('h'); ylabel('Erreur');

% Estimation de l'ordre de convergence
p = polyfit(log(h), log(ErreursL2rel), 1);
ordreL2 = p(1);
loglog(h, exp(p(2))*(h.^p(1)), 'r--');

p = polyfit(log(h), log(ErreursH1rel), 1);
ordreH1 = p(1);
loglog(h, exp(p(2))*(h.^p(1)), 'b--');
fprintf('Ordre de convergence\n\t en norme L2 : %f\n\t en semi-norme H1 : %f\n', ordreL2, ordreH1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2023
