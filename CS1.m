% CS based provenance recovery method on a network of n nodes and h hops
%function []=CS1()
    
clc;clear;
n = 7; % no. of nodes
l = 100; b = 100; %length breadth

% Node distribution
points = [];
for i=1:n
    points = [points; [randi(l,1,1) randi(b,1,1)]];
end

%% Distribution Plot
x = points(:,1);
y = points(:,2);
%figure(1);
%scatter(x,y,7,"blue","filled");
%title('Distrbution of nodes');li
a = [1:n]'; b = num2str(a); c = cellstr(b);
dx = 0.1; dy = 0.1;
text(x+dx, y+dy, c);

%% Distance between coordinates: Adjacency Matrix
Adj_mat = squareform(pdist(points,'euclidean'));
G = graph(Adj_mat);

[T,pred] = minspantree(G); % Form MST using graph G
maxweight = max(T.Edges.Weight); % Max edge weight of MST
reach_mat = Adj_mat; 

% Form connected graph using max edge weight
reach_mat(reach_mat <= maxweight+4 & reach_mat >0) = 1;
reach_mat(reach_mat > maxweight) = 0;
connect_G = graph(reach_mat);

%figure(2);
%plot(connect_G);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation Parameters

E = (n-1)^2;       % Total no of directed edges
EE = n^2;          % Modified for additional constraint
hop_len = [4 5 6]; % Required hop length(to be used later)
h = hop_len(1);    % hop length chosen
%quant_bits = [2 6]; 
mu = 0; sigma = 5;range = 10; % gaussian parameters

error_threshold = 100; % Threshold on error
no_of_pkts = 10^5;     % Total no of pkts txed
error_count = 0;
pkt_count = 0;

mm = 2*h; 
m = [];            % Array for various values of no of rows of matrix A(m,n) 
for ii=1:6
    m = [m mm+(ii-2)*mm/2];
end

% Find optimal path with required hop length 'h'

% Define n nodes with 'class node'
nodes=[]; nodesE=[];  % All varE used for additional constraint
for i= 1:n
    nodes = [nodes node(i,[],reach_mat(i,:))];
    nodesE = [nodesE node(i,[],reach_mat(i,:))];
end

% Destination
dest = nodes(n);
fprintf("Destination node: %d",dest.Node_id);

reach_mat1 = reach_mat;    %adjacency matrix
reach_mat1(reach_mat1 == 0) = inf;
path=[];
[path,src] = find_path(reach_mat1,nodes,h);  %Func to get optimal path like [2,3,1,4]
if isempty(path)
    fprintf('No path available\n\n');
%CS1();
else
    fprintf('Path choosen:');disp(path) 
end

% Path array in terms of 1/0 = edges used/unused
Path_arr = path_array1(path,EE,n);
Path_arrE = path_array1(path,EE,n);
fprintf("path array:");
disp(Path_arr')
%fprintf("path array1:");
%disp(Path_arrE)
%for qBit=1:length(quant_bits)

for m_index =1: length(m)
    error_count = 0;
    error_count_OMP=0;
    %error_rate(qBit,m_index) = 0;
    error_rate(m_index) = 0;
    error_rate_OMP(m_index) = 0;
    pkt_count = 0;
    for i=1:no_of_pkts

        pkt =  Packet;
        %pktE =  Packet;
        Ar = normrnd(mu,sigma,[m(m_index),EE]);

        %ArE = normrnd(mu,sigma,[m(m_index),EE]);
       % A = uencode(Ar,quant_bits(qBit),range,'signed') ; % Quantising Ar

        % Assigning columns to edges
        Edge_id = reshape(Ar,[m(m_index),EE/n,n]);
        %Edge_idE = reshape(ArE,[m(m_index),EE/n,n]);
        pkt_count = pkt_count +1;
        for i= 1:n
            if i==n
                %edge=Ar;
                %edgeE=ArE;
                edge=Edge_id(:,:,i);
            else
                edge=Edge_id(:,:,i);
                %edgeE= Edge_idE(:,:,i);
            end
            nodes(i).Edge_id = edge;
            %nodesE(i).Edge_id = edgeE;
        end

        % Embed edge ids
        pkt = edge_embed_nsquare(m(m_index),path,pkt,nodes);
        %pktE = edge_embed_nsquare(m(m_index),path,pktE,nodesE);
        b = pkt.provenance;
        %bE = pktE.provenance;
        %fprintf("Final provenance\n");disp(b)
        
        % Recovery using OMP
        recovered_x_OMP = OMP(h,b,Ar);

        % Recovery using CVX
        cvx_clear
        cvx_begin quiet
            cvx_precision default
            variable x(EE)
            %x >= 0;
            minimize( norm(x, 1 ) )
            subject to
                %double(A) * x == b;
                Ar * x==b;
                %norm( Ar * x - b) <= e;
        cvx_end


        for k=1:length(x)
            if x(k)<=0.001
                rec_x(k)=0;
            else
               rec_x(k)=1;
            end
        end
        recovered_x = rec_x';
        %fprintf("Path recovered:")
        %disp(recovered_x')
        %newline

        if recovered_x == Path_arr
            %fprintf("Path matched")
            %sprintf('\n');
        else
            error_count = error_count+1; %Increment if path recovered is different from path travelled
        end
        
        if recovered_x_OMP == Path_arr
            %fprintf("Path matched")
            %sprintf('\n');
        else
            error_count_OMP = error_count_OMP +1; %Increment if path recovered is different from path travelled
        end
        
        if error_count == error_threshold % when threshold reached
            %error_rate(qBit,m_index) = error_count/pkt_no;
            %fprintf("Error Rate:%f",error_rate(qBit,m_index));
            break
        end

    end
    %error_rate(qBit,m_index) = error_count/pkt_count;
    error_rate(m_index) = error_count/pkt_count;
    error_rate_OMP(m_index) = error_count_OMP/pkt_count;
    %fprintf("Error Rate:%f for size %d",error_rate(qBit,m_index),quant_bits(qBit)*m(m_index));
    fprintf("Error Rate:%f for column size %d\n",error_rate(m_index),m(m_index));
end
%fprintf("Done for %d bits",quant_bits(qBit));
%prov_size(qBit,:) = quant_bits(qBit)*m; %Provenance size
prov_size = 16*m;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Code for additional constraints
%{ 
%Uncomment only for solving with path constraint
for m_index =1: length(m)
    error_countE = 0;
    %error_rate(qBit,m_index) = 0;
    error_rateE(m_index) = 0;
    pkt_countE = 0;
    for i=1:no_of_pkts

        pktE =  Packet;
        ArE = normrnd(mu,sigma,[m(m_index),EE]);
       % A = uencode(Ar,quant_bits(qBit),range,'signed') ; % Quantising Ar

        % Assigning columns to edges
        Edge_idE = reshape(ArE,[m(m_index),EE/n,n]);
        pkt_countE = pkt_countE +1;
        for i= 1:n
           edgeE= Edge_idE(:,:,i);
           nodesE(i).Edge_id = edgeE;
        end

        % Embed edge ids

        pktE = edge_embed_nsquare(m(m_index),path,pktE,nodesE);
        bE = pktE.provenance;
        %fprintf("Final provenance\n");disp(b)

        % Recovery
        % With n^2 values in x and additional constraint

        B=[];
        for ii=1:n
            count=0;
            for j=1:EE
                k=n*(ii-1)+1;
                count=k+n;
                if j>=k && j<count 
                    B(ii,j)=1;
                else 
                    B(ii,j)=0;
                end
            end
        end
        C=[];
        for ii=1:n
            C=[C;eye(n)];
        end

        cvx_begin quiet
            cvx_precision default
            variable xE(EE)
            %A_mat = xE*xE';
            %Z = B*A_mat*C;
            %Zt = Z^h;
            %flag = any(Zt(end,:));
            minimize( norm(xE, 1 ) )
            subject to
                ArE * xE == bE;
             %   flag == 1;
        cvx_end

        for k=1:length(xE)
            if xE(k)<=0.001
                rec_xE(k)=0;
            else
               rec_xE(k)=1;
            end
        end
        recovered_xE = rec_xE';

        A_mat = rec_xE*rec_xE';
        Z = B*A_mat*C;
        Zt = Z^h;
        flag = any(Zt(end,:));
        if flag==1
            %fprintf("X is:\n")
            %fprintf("Valid")
            %disp(xE)
        else
            pkt_countE = pkt_countE -1;
            continue
        end

        if recovered_xE == Path_arrE
            %fprintf("Path matched")
            sprintf('\n');
        else
            error_countE = error_countE+1; %Increment if path recovered is different from path travelled
        end
        if error_countE == error_threshold % when threshold reached
            %error_rate(qBit,m_index) = error_count/pkt_no;
            %fprintf("Error Rate:%f",error_rate(qBit,m_index));
            break
        end
    end
    %error_rate(qBit,m_index) = error_count/pkt_count;
    error_rateE(m_index) = error_countE/pkt_countE;
    %fprintf("Error Rate:%f for size %d",error_rate(qBit,m_index),quant_bits(qBit)*m(m_index));
    fprintf("Error Rate w constraint:%f for column size %d\n",error_rateE(m_index),m(m_index));
end

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%plot(prov_size,error_rate);
figure(1)
semilogy(m,error_rate, 'mo-', 'LineWidth', 2);  % With cvx
hold on
%semilogy(m,error_rateE, 'b+-', 'LineWidth', 2);
semilogy(m,error_rate_OMP, 'b+-', 'LineWidth', 2);  % with OMP
axis([0 120 0 1]);
grid on
legend('CVX', 'OMP');
title('Error rate vs provenance size');
xlabel('Column size');
ylabel('Error rate')
% hold off
%end
