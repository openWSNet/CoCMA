%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [coverage_rec,avg_packets_to_bs,avg_packets_to_ch,DEAD,S,last_round,CLUSTERHS,avg_ch]=LEACH(sense_node,...
                                              rmax, p, best_gene, node_x, node_y, sink_x, sink_y,...
                                              packet_bit, local_search_for_alive_neighbor_flag) 
%avg_ch: average of culser heads when all nodes run out
global target_x target_y grid_range_x grid_range_y span sense_range;

if local_search_for_alive_neighbor_flag==1
    gene_of_each_round=best_gene;
else
    gene_of_each_round=ones(1,length(best_gene));
end

Eo=0.25;
%Eelec=Etx=Erx
ETX=50*0.000000001;
ERX=50*0.000000001;
%Transmit Amplifier types
Emp=0.1*0.000000001;
%Data Aggregation Energy
EDA=5*0.000000001;
flag=0;                                                         
%Data Aggregation Energy

%Computation of do

for i=1:1:sense_node
        S(i).xd=node_x(i);
        XR(i)=S(i).xd;
        S(i).yd=node_y(i);
        YR(i)=S(i).yd;
        S(i).G=0;
        S(i).sons=0;
        %initially there are no cluster heads only nodes
        S(i).type='N';
        S(i).E=Eo;    
end

S(sense_node+1).xd=sink_x;
S(sense_node+1).yd=sink_y;

dist_node=zeros(sense_node+1,sense_node+1); % calculate the distance between two nodes¡Asense_node+1 indicates of sink
for i=1:(sense_node+1)
    for j=1:(sense_node+1)
        dist_node(i,j)=dist(S(i).xd,S(i).yd,S(j).xd,S(j).yd);
    end
end

coverage_rec=zeros(1,rmax+1);


countCHs=0;
%counter for CHs per round
rcountCHs=0;
cluster=1;

countCHs;
rcountCHs=rcountCHs+countCHs;
flag_first_dead=0;
last_round=0;

global_count=0;
%%%%%%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%% START OF ROUNDS %%%%%%%%%%%%%%%%%%%%%%%%%%
for r=1:1:rmax
  %Operation for epoch
  r
    if mod(r, round(1/p) )==0
       for i=1:sense_node  % update G, in order to determine if the epoch is done
         S(i).G=0;
         S(i).cl=0;
       end
    end
    for y=1:sense_node
       S(y).sons=0; % record the number of nodes belonging to each head
    end
% hold off;

%Number of dead nodes
dead=0;
%Number of dead Advanced Nodes
dead_a=0;
%Number of dead Normal Nodes
dead_n=0;

%counter for bit transmitted to Bases Station and to Cluster Heads
packets_TO_BS=0;
packets_TO_CH=0;
%counter for bit transmitted to Bases Station and to Cluster Heads 
%per round
PACKETS_TO_CH(r+1)=0;
PACKETS_TO_BS(r+1)=0;


for i=1:1:sense_node
    %checking if there is a dead node
    if (S(i).E<=0)
%         disp('dead=');
        dead=dead+1;
        if local_search_for_alive_neighbor_flag==1  && gene_of_each_round(i)==1
            global_count=global_count+1;
%             clf;
            gene_of_each_round_old=gene_of_each_round; % record the old topology
            gene_of_each_round(i)=0;  % update the gene and re-evaluate if the node is dead                                
            gene_of_each_round=local_search_for_alive_neighbor(i,gene_of_each_round,S,dist_node);
%             plot_topology(gene_of_each_round_old,gene_of_each_round,sense_range,node_x,node_y,target_x,target_y);
%             pause;         
        end 
            gene_of_each_round(i)=0;  %double check 
            S(i).type='D';          
    end
  
end

% recheck 
if sum(gene_of_each_round)==0
    for b=1:sense_node
        if S(b).E>0
            gene_of_each_round(b)=1;
        end
    end
end

for j=1:sense_node
    if S(j).E>0 && gene_of_each_round(j)==1
        S(j).type='N';
    end
end

% after check, evaluate the coverage
 [coverage_rec(1,r+1),temp]=fit_foreach(gene_of_each_round);

STATISTICS(r+1).DEAD=dead;
DEAD(r+1)=dead;


%When the first node dies
if (dead==1)
    if(flag_first_dead==0)
        first_dead=r
        flag_first_dead=1;
    end
end

countCHs=0; 
cluster=1; %counter index
alive_node=0;
for i=1:1:sense_node
   if(S(i).E>0)
   temp_rand=rand;    
   alive_node=alive_node+1;
   if ( (S(i).G)<=0  && gene_of_each_round(i)==1) % possibly act as a head if G<=0
       %Election of Cluster Heads    
       if(temp_rand<= (p/(1-p*mod(r,round(1/p)))))
            countCHs=countCHs+1;
            packets_TO_BS=packets_TO_BS+1;  
            S(i).type='C';
            S(i).G=round(1/p)-1;   
            % G will not be 0 if it has acted as a head. G can be used as a index to determine if act as a head
            C(cluster).xd=S(i).xd;
            C(cluster).yd=S(i).yd;
            
            distance=dist_node(i,sense_node+1);% distance between the node and sink
            C(cluster).distance=distance;
            C(cluster).id=i;
            X(cluster)=S(i).xd;
            Y(cluster)=S(i).yd;
            cluster=cluster+1;
        end     
    
    end
   end 

end

CLUSTERHS(r+1)=(cluster-1);
if (alive_node==0 && flag==0 )||(r==rmax && flag==0 ) % record the last_round if the max round is achieved or all nodes are dead
    last_round=r;
    flag=1;
end
packets_TO_CH=0;
%Election of Associated Cluster Head for Normal Nodes
for i=1:1:sense_node
   if ( S(i).type=='N' && S(i).E>0 && gene_of_each_round(i)==1)
     if(cluster-1>=1) % if there are more than one head
       min_dis=dist_node(i,sense_node+1);
       min_dis_cluster=1;
       for c=1:1:cluster-1
           temp=min(min_dis,sqrt( (S(i).xd-C(c).xd)^2 + (S(i).yd-C(c).yd)^2 ) );
           if ( temp<min_dis )
               min_dis=temp;
               min_dis_cluster=c;
           end
       end     
          S(i).E=S(i).E- ( ETX*(packet_bit) + Emp*packet_bit*( min_dis * min_dis)); 
        if(min_dis_cluster~=1)
          S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E- ( (ERX + EDA)*packet_bit ); 
          S(C(min_dis_cluster).id).sons=S(C(min_dis_cluster).id).sons+1;
          packets_TO_CH=packets_TO_CH+1;
        end
        S(i).min_dis=min_dis;
        S(i).min_dis_cluster=min_dis_cluster;      
     else
        min_dis=dist_node(i,sense_node+1);
        min_dis_cluster=1;
        S(i).E=S(i).E- ( ETX*(packet_bit) + Emp*packet_bit*( min_dis * min_dis)); 
        packets_TO_BS=packets_TO_BS+1; 
     end
 end
end
      PACKETS_TO_CH(r+1)=packets_TO_CH; 
      PACKETS_TO_BS(r+1)=packets_TO_BS;     
         % the number of packets to be sent to BS from heads
for h=1:sense_node %% calculate the energy consumption of the transmission between head and bs
    if S(h).type=='C' && S(h).E>0
       dis=dist_node(h,sense_node+1);% distance between the node and sink
       S(h).E=S(h).E- ( EDA*(packet_bit) + ETX*(packet_bit) + Emp*packet_bit*( dis * dis)); 
    end
end
% hold on;

countCHs;
rcountCHs=rcountCHs+countCHs;

for i=1:sense_node %% reset the energy index to 0 if it's minus 
    if S(i).E<0
        S(i).E=0;
    end
end

% if r==1 || r==500 || r==1500 || r==3000 || r==4500 || r==4900
% % if r==rmax
%     energy_analysis(100,100,S,sense_node);  
% %     print(figure(10),'-dmeta',[int2str(r)]);  
% 
%     pause;
% end

end

% energy_analysis(grid_range_x*span,grid_range_y*span,S,sense_node);
% plot_topology(1,sense_range,gene_of_each_round,node_x,node_y,target_x,target_y);
gene_of_each_round
fprintf('\nglobal_count=%d',global_count);

figure(2);
if last_round~=0
    avg_ch=sum(CLUSTERHS)/last_round;
    avg_packets_to_bs=sum(PACKETS_TO_BS)/last_round;
    avg_packets_to_ch=sum(PACKETS_TO_CH)/last_round;
else
    avg_ch=sum(CLUSTERHS)/rmax;
end

plot(sense_node-DEAD(2:(rmax+1)),'b');     %  the dead nodes will appear from 2nd round 
figure(3);
target_num=length(target_x(1,:))*length(target_y(:,1));

plot(coverage_rec(2:(rmax+1))/target_num,'r')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   STATISTICS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                     %
%  DEAD  : a rmax x 1 array of number of dead nodes/round 
%  DEAD_A : a rmax x 1 array of number of dead Advanced nodes/round
%  DEAD_N : a rmax x 1 array of number of dead Normal nodes/round
%  CLUSTERHS : a rmax x 1 array of number of Cluster Heads/round
%  PACKETS_TO_BS : a rmax x 1 array of number packets send to Base Station/round
%  PACKETS_TO_CH : a rmax x 1 array of number of packets send to ClusterHeads/round
%  first_dead: the round where the first node died                   
%                                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end




