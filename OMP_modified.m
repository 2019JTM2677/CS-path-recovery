% Function for modified OMP algorithm 
% K sparsity, x recovered signal, y measured signal, A matrix
function [x] = OMP_modified (K,y,A)

[m,n] = size (A) ;
x = zeros (n,1) ;
L=8;    % No. of columns chosen in each iteration
Residue = zeros(m,L,K);
kk = zeros(L,K);
ind = [];
val = [];
C = [];V=[];

% Iterating K times 
    for J = 1 : K
        x_cap=[];
        %Index Search
        if J==1
            [~ ,index] = maxk(abs(A' * y),L); %choose first L max elemnts
            %index = sort(index);
            %v = sort(v);
            %display(abs(A' * y));
            %display(maxk(abs(A' * y),5));
        else
            for li=1:L
                [val(:,li),ind(:,li)] = maxk(abs(A' * Residue(:,li,J-1)),L);
                %display(maxk(abs(A' * Residue(:,li,J-1)),L));
                %display(maxk(abs(A' * y),5));
            end
            ind1 = reshape(ind,[L^2,1]);
            val1 = reshape(val,[L^2,1]);
            [~,maxind] = maxk(val1,L);
            maxind = sort(maxind);
            V=C;
            kv = kk;
            
            for li=1:L
                r = ceil(maxind(li)/L);
                if r==li
                    %fprintf("leave previous basis matrix same for that li")
                else  %Change the basis matrix
                    V(:,:,li) = C(:,:,r);
                    kv(li,:) = kk(r,:) ;
                end
            end
            C = V;
            kk = kv;
            
            
            index = ind1 (maxind);
        end
            kk (:,J) = index ;  % support for L possible x
            l = index;          %list of columns selected in each iteration
        
        % For each element in l find x_cap 
        B=[];
        for li = 1:L
                      
            B = cat(3,B,A(:,l(li)));  %columns to be added for different L
            
            %fprintf("matrix is:")
            %disp(B'*B);
            %fprintf("\n");
            %B_mod = inv(B'*B);
                        
        end
        C = cat(2,C,B);   % Basis matrix
        %display(C);
        for li = 1:L
            %fprintf("li=%d",li);
            C_mod = C(:,:,li)'*C(:,:,li);
            y_mod = C(:,:,li)'*y;
            x_cap(:,li) = C_mod\y_mod;%B'*y  ;         % x^ in each iteration
            %fprintf("x cap= \n");
            %display(C(:,:,li));
            %display(x_cap(:,li));
            %fprintf("x=%d\n",x_cap(:,li));
            Residue(:,li,J) = y - (C(:,:,li)*(x_cap(:,li)));%new_x';  % residue 
        end
    end
    
    for li=1:L
        for i=1:size(kk,2)
            x(kk(li,i),li)= x_cap(i,li);     % Final x after recovery
        end
    end
end