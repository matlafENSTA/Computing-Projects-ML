function u = householder_vector(x)
    u = x;
    u(1) = u(1) + sign(u(1)) * norm(u);
end
