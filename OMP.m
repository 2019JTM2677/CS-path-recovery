% Function for OMP algorithm 
% K sparsity, x recovered signal, y measured signal, A matrix
function [x] = OMP (K,y,A)

[m,n] = size (A) ;
x = zeros (n,1) ;
Residue = zeros(m,K);
B=[];
x_cap=zeros(K,1);
kk = zeros();

% Iterating K times 
for J = 1 : K
    
    %Index Search
    if J==1
        [~ ,index] = max(abs(A' * y)) ;
    else
        [~,index] = max(abs(A' * Residue(:,J-1)));
    end
        kk (J) = index ;
    
    %Residue Update
    %w (:,J) = A (:,kk (J)) ;    % Choosing column of A for basis matrix
    B = [B A(:,kk(J))];          % Basis matrix
    B_mod = inv(B'*B);
    x_cap = B_mod*B'*y  ;         % x^ in each iteration
    Residue(:,J) = y - B*x_cap;  % residue 
    
end
for i=1:size(kk,2)
    x(kk(i))= x_cap(i);     % Final x after recovery
end
