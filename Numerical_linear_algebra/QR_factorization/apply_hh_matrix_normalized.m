function A = apply_hh_matrix_normalized(A, u)
    % un bon point de départ est le code pour apply_hh_matrix!
    [m, n] = size(A);
    assert(length(u) == m - 1, 'La longueur du vecteur u doit être égale à m - 1');
    beta = 2 / (u * u' + 1);
    for j = 1:n
        col = A(:, j);
        coeff = [1, u'] * col * beta;
        col(1) = col(1) - coeff;
        for i = 1:m-1
            col(i) = col(i) - coeff * u(i);
        end
        A(:, j) = col;
    end
end
