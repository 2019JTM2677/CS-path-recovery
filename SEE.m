% Single edge embedding on a network of n nodes and h hops
%function []= SEE()
    clc;clear;
    n = 5; % no. of nodes
    hop_len = [4 5 6]; % Required hop length
    h= hop_len(1); % choose hop length according to path
    % Define path to be travelled
    path = [1 3 2 4 n];
    
    
    reach_mat = ones(n) - diag(ones([1,n])); %adjacency matrix for complete graph
    connect_G = graph(reach_mat~=0);
    %figure(2);
    %plot(connect_G);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Simulation
    % Simulation Parameters
    E = (n-1)^2; % Total no of directed edges
    mu=0;sigma=5; %gaussian parameters
    error_threshold = 100; % Error Rate = error_threshold/no_of_pkts
    no_of_pkts = 10^5; % Total no of pkts txed
    
    
    mm = 2*h; 
    m=[]; % form array for various values of no of rows of matrix A(m,n) 
    for ii=1:7
        m=[m mm+(ii-2)*mm/2];
    end
    
    % Define n nodes 
    nodes=[]; 
    for i= 1:n
        nodes=[nodes node(i,[],reach_mat(i,:))];
    end
    
    % Destination
    dest = nodes(n);
    %fprintf("Destination node: %d",dest.Node_id);
    
    src = nodes(path(1));
    fprintf('Path choosen:');disp(path) 
    
    % Path array similar to x (y=Ax) for verify
    Path_arr = path_array(path,E,n);
    %fprintf("\npath array:");
    %disp(Path_arr)
    
    fprintf("\n Single edge embedding\n");
    
    for m_index =1: length(m)
        error_count = 0;
        error_count_OMP=0;
        error_rate(m_index) = 0;
        error_rate_OMP(m_index) = 0;
        pkt_count = 0;

        for i=1:no_of_pkts

            pkt =  Packet;
            pkt_count = pkt_count +1;
            Ar = normrnd(mu,sigma,[m(m_index),E]);

            % Assigning columns to edges
            Edge_id = reshape(Ar,[m(m_index),E/(n-1),n-1]);

            for ni= 1:n
                if ni==n
                    edge=Ar;
                else
                    edge=Edge_id(:,:,ni);
                end
                nodes(ni).Edge_id = edge;
            end

            % Embed edge ids
            pkt = edge_embed(m(m_index),path,pkt,nodes);
            b = pkt.provenance;
            %fprintf("Final provenance\n");disp(b)

            % Recovery using OMP
            x_OMP = OMP(h,b,Ar);
            for k=1:length(x_OMP)
                if abs(x_OMP(k))<=0.001
                    rec_x_OMP(k)=0;
                else
                   rec_x_OMP(k)=1;
                end
            end
            
            if rec_x_OMP' == Path_arr
            %    fprintf("Path matched\n")
            else
               error_count_OMP = error_count_OMP +1; % Increment count
            end

            % Recovery using CVX
            x = cvx_solver(E,b,Ar);
            for k=1:length(x)
                if abs(x(k))<=0.001
                    rec_x(k)=0;
                else
                   rec_x(k)=1;
                end
            end
            % Comparing recovered path wth original path
            if rec_x' == Path_arr
                %fprintf("Path matched\n")
            else
                error_count = error_count+1; %Increment if path recovered is different from path travelled
            end
            if error_count == error_threshold % when threshold reached
                break
            end
        end

        error_rate(m_index) = error_count/pkt_count;
        error_rate_OMP(m_index) = error_count_OMP/pkt_count;
        fprintf("Error Rate:%f for column size %d\n",error_rate(m_index),m(m_index));
    end
        
    
    
    figure(1)
    semilogy(m,error_rate, 'mo-', 'LineWidth', 2);
    hold on
    semilogy(m,error_rate_OMP, 'bo-', 'LineWidth', 2);
    axis([0 120 0 1]);
    grid on
    %legend('SE CVX', 'SE OMP');
    title('Error rate vs number of rows');
    xlabel('Column size');
    ylabel('Error rate')
    %hold off
%end
DEE(path,n)