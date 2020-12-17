% CS based provenance recovery method on a network of n nodes and h hops
% Comparing SE and DE
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation Parameters
%E = (n-1)^2; % Total no of directed edges
EE = n^2;    % Modified for additional constraint
DE = ((n-2)^2)*(n-1); % No. of Double edges  
hop_len = [3 4 5 6]; % Required hop length(to be used later)
h = hop_len(2);    % hop length chosen
mu = 0; sigma = 5; %gaussian parameters
range = 12;
error_threshold = 300; % Error rate= error_threshold/ no_of_pkts
no_of_pkts = 1000 ;% Total no of pkts txed

mm = 2*h; 
m = []; % form array for various values of no of rows of matrix A(m,n) 
for ii=1:7
    m = [m mm+(ii-2)*mm/2];
end

% Finding optimal path with required hop length

% Define n nodes 
nodes=[]; nodes_de=[];  % All varE used for additional constraint
for i= 1:n
    nodes = [nodes node(i,[],reach_mat(i,:))];
    nodes_de = [nodes_de node(i,[],reach_mat(i,:))];
end

% Destination
dest = nodes(n);
fprintf("Destination node: %d",dest.Node_id);

reach_mat1= reach_mat;
reach_mat1(reach_mat1 == 0) = inf;
path=[];
[path,src] = find_path(reach_mat1,nodes,h);
if isempty(path)
fprintf('No path available\n\n');
%CS1();
else
fprintf('\nPath choosen:');disp(path) 
end

% Path array similar to x (y=Ax) for verify
Path_arr = path_array1(path,EE,n);
Path_arr_de = path_array_de(path,DE,n);
%fprintf("path array:");
%disp(Path_arr')
fprintf("path array1:");
disp(Path_arr_de')

for m_index =1: length(m)
    error_count = 0;
    error_count_OMP=0;
    error_count_de=0;
    error_count_Ode=0;
    
    error_rate(m_index) = 0;
    error_rate_OMP(m_index) = 0;
    error_rate_de(m_index) = 0;
    error_rate_Ode(m_index) = 0;
    
    pkt_count = 0;
    
    for i=1:no_of_pkts

        pkt =  Packet;
        pkt_de =  Packet;
        Ar = normrnd(mu,sigma,[m(m_index),EE]);
        Ar_de = normrnd(mu,sigma,[m(m_index),DE]);
        %ArE = normrnd(mu,sigma,[m(m_index),EE]);
       % A = uencode(Ar,quant_bits(qBit),range,'signed') ; % Quantising Ar

        % Assigning columns to edges
        Edge_id = reshape(Ar,[m(m_index),EE/n,n]);
        Edge_id_de = reshape(Ar_de,[m(m_index),(n-2)^2,n-1]);
        pkt_count = pkt_count +1;
        
        for ni= 1:n
            if ni==n
                %edge=Ar;
                edge = Edge_id(:,:,ni);
                dbl_edge = Ar_de;
            else
                edge = Edge_id(:,:,ni);
                dbl_edge = Edge_id_de(:,:,ni);
            end
            nodes(ni).Edge_id = edge;
            nodes_de(ni).Edge_id = dbl_edge;
        end

        % Embed edge ids
        pkt = edge_embed_nsquare(m(m_index),path,pkt,nodes);
        pkt_de = double_edge_embed(m(m_index),path,pkt_de,nodes_de);
        b = pkt.provenance;
        b_de = pkt_de.provenance;
        %fprintf("Final provenance\n");disp(b_de)
        
        % Recovery using OMP
        x_OMP = OMP(h,b,Ar);
        rec_x_OMP = uint8(x_OMP);
        if rec_x_OMP == Path_arr
        %    fprintf("Path matched\n")
        else
           error_count_OMP = error_count_OMP +1; % Increment count
        end
        
       % Recovery using OMP double edge
       if rem(length(path),2)==0
           h_omp = length(path)/2 -1;
       else
           h_omp = floor(length(path)/2);
       end
       x_Ode = OMP(h_omp,b_de,Ar_de);
       rec_x_Ode = uint8(x_Ode);
       if rec_x_Ode == Path_arr_de
           %flag
       else
           error_count_Ode = error_count_Ode +1; % Increment count
          %  flag_OMP=0;
       end

       
        % Recovery using CVX
        x = cvx_solver(EE,b,Ar);
        rec_x = uint8(x);
        %fprintf("Path recovered:")
        %disp(recovered_x')
        %newline

        if rec_x == Path_arr
            %fprintf("Path matched")
            %sprintf('\n');
        else
            error_count = error_count+1; %Increment if path recovered is different from path travelled
        end
        
        % Recovery using CVX double edge
        
        x_de = cvx_solver(DE,b_de,Ar_de);
        rec_x_de = uint8(x_de);
        
        if rec_x_de == Path_arr_de
            %fprintf("Path matched")
            %sprintf('\n');
        else
            error_count_de = error_count_de +1; %Increment if path recovered is different from path travelled
        end

        
        if error_count == error_threshold% when threshold reached
            break
        end

    end
    %error_rate(qBit,m_index) = error_count/pkt_count;
    error_rate(m_index) = error_count/pkt_count;
    error_rate_OMP(m_index) = error_count_OMP/pkt_count;
    error_rate_de(m_index) = error_count_de/pkt_count;
    error_rate_Ode(m_index) = error_count_Ode/pkt_count;
    fprintf("Error Rate:%f for column size %d\n",error_rate(m_index),m(m_index));
end
% end

%%
%plot(prov_size,error_rate);
figure(1)
a1=semilogy(m,error_rate, 'm*-', 'LineWidth', 2);
hold on
a2=semilogy(m,error_rate_OMP, 'bs-', 'LineWidth', 2);
a3=semilogy(m,error_rate_de, 'g*-', 'LineWidth', 2);
a4=semilogy(m,error_rate_Ode, 'rs-', 'LineWidth', 2);
axis([0 120 0 1]);
grid on
%M1="n=5,h=4"; M2="n=9,h=4";
legend('SE CVX', 'SE OMP','DE CVX', 'DE OMP');
%legend([a1,a2], [M1, M2]);
title('Error rate vs provenance size');
xlabel('Column size, m');
ylabel('Reconstruction Error rate');
hold off
%end