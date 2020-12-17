% Function for CVX solver 
% E dimension of x, x recovered signal, y measured signal, A matrix
function [x] = cvx_solver (E,y,A)
    cvx_begin quiet
            cvx_precision default
            variable x(E)
            minimize( norm(x, 1 ) )
            subject to
                A*x == y;
        cvx_end

end
