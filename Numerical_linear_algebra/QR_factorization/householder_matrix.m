function F = householder_matrix(u)
    F = eye(length(u)) - 2 * (u * u') / (u' * u);
end
