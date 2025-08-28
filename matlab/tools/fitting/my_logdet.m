function logdet = my_logdet(A)
% calculate the log determinant of a matrix
% A: a square matrix

% [L,U] = lu(A);
% logdet = sum(log(diag(U)));
% logdet = sum(log(abs(diag(U))));
% tmp = diag(L).*diag(U);


% logdet = trace(log(U));

% return

% cholesky decomposition
% suppose the fastest way for a positive definite matrix
logdet = 2 * sum(log(diag(chol(A))));
% e = eig(A);
end