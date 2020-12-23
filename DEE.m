% Double edge embedding
function []= DEE(path,n)
    %n = 13; % no. of nodes

    reach_mat = ones(n) - diag(ones([1,n]));
    connect_G = graph(reach_mat~=0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Simulation Parameters
    DE = ((n-2)^2)*(n-1); % No. of Double edges  
    hop_len = [4 5 6]; % Required hop length(to be used later)
    h = hop_len(1);     %  hop length chosen
    mu = 0; sigma = 5;  % gaussian parameters

    error_threshold = 100; % Error rate= error_threshold/ no_of_pkts
    no_of_pkts = 1000 ;% Total no of pkts txed

    mm = 2*h; 
    m = []; % form array for various values of no of rows of matrix A(m,n) 
    for ii=1:7
        m = [m mm+(ii-2)*mm/2];
    end
    
    % Define n nodes
    nodes_de=[];  % All varE used for additional constraint
    for i= 1:n
        nodes_de = [nodes_de node(i,[],reach_mat(i,:))];
    end

    % Destination
    dest = nodes_de(n);
    %fprintf("Destination node: %d",dest.Node_id);

    %Define path
    %path = [1 3 2 4 13];
    src = nodes_de(path(1));
    %fprintf('Path choosen:');disp(path) 
    
    % Path array similar to x (y=Ax) for verify
    Path_arr_de = path_array_de(path,DE,n);
    %fprintf("path array for double edge:");
    %disp(Path_arr_de')

    fprintf("\n Double edge embedding\n");
    for m_index =1: length(m)
        error_count_de=0;
        error_count_Ode=0;

        error_rate_de(m_index) = 0;
        error_rate_Ode(m_index) = 0;

        pkt_count = 0;

        for i=1:no_of_pkts
            pkt_de =  Packet;
            Ar_de = normrnd(mu,sigma,[m(m_index),DE]);
            
            % Assigning columns to edges
            Edge_id_de = reshape(Ar_de,[m(m_index),(n-2)^2,n-1]);
            pkt_count = pkt_count +1;

            for ni= 1:n
                if ni==n
                    dbl_edge = Ar_de;
                else
                    dbl_edge = Edge_id_de(:,:,ni);
                end
                nodes_de(ni).Edge_id = dbl_edge;
            end

            % Embed double edge ids
            pkt_de = double_edge_embed(m(m_index),path,pkt_de,nodes_de);
            b_de = pkt_de.provenance;
            %fprintf("Final provenance\n");disp(b_de)

           % Recovery using OMP double edge
           if rem(length(path),2)==0    
               h_omp = length(path)/2 -1;
           else
               h_omp = floor(length(path)/2);
           end
           
           x_Ode = OMP(h_omp,b_de,Ar_de);
           rec_x_Ode = uint8(x_Ode);
           
           if rec_x_Ode == Path_arr_de
               %fprintf("Path matched")
           else
               error_count_Ode = error_count_Ode +1; % Increment count
           end


            % Recovery using CVX double edge
            x_de = cvx_solver(DE,b_de,Ar_de);
            rec_x_de = uint8(x_de);

            if rec_x_de == Path_arr_de
                %fprintf("Path matched\n")
            else
                error_count_de = error_count_de +1; %Increment if path recovered is different from path travelled
            end
            
            if error_count_de == error_threshold% when threshold reached
                break
            end

        end
        error_rate_de(m_index) = error_count_de/pkt_count;
        error_rate_Ode(m_index) = error_count_Ode/pkt_count;
        fprintf("Error Rate:%f for column size %d\n",error_rate_de(m_index),m(m_index));
    end
    % end

    %% PLOTS
    figure(1)
    %a1=semilogy(m,error_rate, 'm+-', 'LineWidth', 2);
    %a2=semilogy(m,error_rate_OMP, 'bd-', 'LineWidth', 2);
    a3=semilogy(m,error_rate_de, 'gd-', 'LineWidth', 2);
    hold on
    %a4=semilogy(m,error_rate_Ode, 'ro-', 'LineWidth', 2);
    axis([0 120 0 1]);
    grid on
    %M1="n=5,h=4"; M2="n=9,h=4";
    legend('SE CVX','DE CVX');
    %legend([a1,a2], [M1, M2]);
    title('Error rate vs provenance size');
    xlabel('Column size, m');
    ylabel('Reconstruction Error rate');
    %hold off
end