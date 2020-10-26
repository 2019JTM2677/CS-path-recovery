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
    x=points(:,1);
    y=points(:,2);
    %figure(1);
    %scatter(x,y,7,"blue","filled");
    %title('Distrbution of nodes');li
    a = [1:n]'; b = num2str(a); c = cellstr(b);
    dx = 0.1; dy = 0.1;
    text(x+dx, y+dy, c);
    
    %% Distance between coordinates: Adjacency Matrix
    Adj_mat = squareform(pdist(points,'euclidean'));
    G = graph(Adj_mat);
    
    %plot(G,'EdgeLabel',G.Edges.Weight);
    %p=plot(G);
    %highlight(p,T)
    %plot(T,'EdgeLabel',T.Edges.Weight);
    
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
    %% Simulation
    % Simulation Parameters
    E = (n-1)^2; % Total no of directed edges
    hop_len = [5 4 6]; % Required hop length
    %quant_bits = [2 6]; 
    mu=0;sigma=5; %gaussian parameters
    range=10;
    error_threshold = 10; % Threshold on error
    no_of_pkts = 10^3; % Total no of pkts txed
    error_count = 0;
    pkt_count = 0;
    
    mm=2*hop_len(2); 
    m=[]; % form array for various values of no of rows of matrix A(m,n) 
    for ii=1:6
        m=[m mm+(ii-2)*mm/2];
    end
    
    % Finding optimal path with required hop length
    
    % Define n nodes 
    nodes=[]; 
    for i= 1:n
        nodes=[nodes node(i,[],reach_mat(i,:))];
    end
    
    % Destination
    dest = nodes(n);
    fprintf("Destination node: %d",dest.Node_id);
    
    reach_mat1= reach_mat;
    reach_mat1(reach_mat1 == 0) = inf;
    path=[];
    [path,src] = find_path(reach_mat1,nodes,hop_len(2));
    if isempty(path)
        fprintf('No path available\n\n');
        %CS1();
    else
        fprintf('Path choosen:');disp(path) 
    end
    
    % Path array similar to x (y=Ax) for verify
    Path_arr = path_array(path,E,n);
    fprintf("path array:");
    disp(Path_arr)
    
    %for qBit=1:length(quant_bits)
        for m_index =1: length(m)
            error_count = 0;
            %error_rate(qBit,m_index) = 0;
            error_rate(m_index) = 0;
            pkt_count = 0;
            for i=1:no_of_pkts
                
                pkt =  Packet;
                pkt_count = pkt_count +1;
                Ar = normrnd(mu,sigma,[m(m_index),E]);
                
               % A = uencode(Ar,quant_bits(qBit),range,'signed') ; % Quantising Ar
                
                % Assigning columns to edges
                Edge_id = reshape(Ar,[m(m_index),E/(n-1),n-1]);
                
                for i= 1:n
                    if i==n
                        edge=Ar;
                    else
                        edge=Edge_id(:,:,i);
                    end
                    nodes(i).Edge_id = edge;
                end
                
                % Embed edge ids
                pkt = edge_embed(m(m_index),path,pkt,nodes);
                b = pkt.provenance;
                %fprintf("Final provenance\n");disp(b)
                
                % Recovery
                %e =  0.001;
                cvx_begin quiet
                    cvx_precision default
                    variable x(E)
                    %x >= 0;
                    minimize( norm(x, 1 ) )
                    subject to
                        Ar * x==b;
                        %norm( Ar * x - b) <= e;
                cvx_end
                 
                % x found has -0.000 as elements, for comparison changing that to 0
                for k=1:length(x)
                    if x(k)<=0.001
                        rec_x(k)=0;
                    else
                       rec_x(k)=1;
                    end
                end
                recovered_x = rec_x';
                
                % Comparing recovered path wth original path
                if recovered_x == Path_arr
                    %fprintf("Path matched")
                    sprintf('\n');
                else
                    error_count = error_count+1; %Increment if path recovered is different from path travelled
                end
                if error_count == error_threshold % when threshold reached
                    %error_rate(qBit,m_index) = error_count/pkt_no;
                    %fprintf("Error Rate:%f",error_rate(qBit,m_index));
                    break
                end
            end
            %error_rate(qBit,m_index) = error_count/pkt_count;
            error_rate(m_index) = error_count/pkt_count;
            %fprintf("Error Rate:%f for size %d",error_rate(qBit,m_index),quant_bits(qBit)*m(m_index));
            fprintf("Error Rate:%f for column size %d",error_rate(m_index),m(m_index));
        end
        
        prov_size = 16*m;
   % end
    
    
    figure(1)
    semilogy(m,error_rate, 'mo-', 'LineWidth', 2);
    axis([0 120 0 1]);
    grid on
    title('Error rate vs number of rows');
    xlabel('Column size');
    ylabel('Error rate')
    
%end