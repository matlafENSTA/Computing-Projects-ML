FUNCTION [Q, R] = HOUSEHOLDER_QR(A)
    R = A;
    [M, N] = SIZE(A);
    KMAX = MIN(M, N);
    QADJ = EYE(M, M);  % MATRICE IDENTITÉ DE TAILLE M×M
    FOR K = 1:KMAX
        Z = ZEROS(M,M);
        X = R(K:M, K);
        U = HOUSEHOLDER_VECTOR(X);
        F = HOUSEHOLDER_MATRIX(U)
        Z(M-KMAX+K:M,M-KMAX+K:M) = F;
        QADJ = QADJ * Z
        R(K:M, :) = F * R(K:M, :)
    END

    Q = CONJ(QADJ');  % ADJOINT(QADJ) EN JULIA EST CONJ(QADJ') EN MATLAB
END
