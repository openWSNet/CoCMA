%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DEAD,avg_packets_to_bs,avg_packets_to_ch,last_round,avg_ch]=LEACH_original_m(rmax,rs,p) 
% load data1;     %load the same node_x node_y

sense_node=400; 
packet_bit=2000;
grid_range_x=200;
grid_range_y=200;

span=0.04;      
sink_x=50;
sink_y=200;
target_x=zeros(grid_range_y*span,grid_range_x*span);
target_y=zeros(grid_range_y*span,grid_range_x*span);
dist_node_target=zeros(grid_range_y*span,grid_range_x*span,sense_node);
% for k=1:sense_node
%        node_x(k)=fix(rand*grid_range_x*span);
%        node_y(k)=fix(rand*grid_range_y*span);
% end
m=0;
n=0;
for k=1:400
    node_x(k)=m;
    node_y(k)=n;
    if m>=95
        m=0;
        n=n+5;
    else
        m=m+5;
    end
end

for i=1:grid_range_y*span  
    for j=1:grid_range_x*span
        target_x(i,j)=6.25+(j-1)*12.5;
        target_y(i,j)=6.25+(i-1)*12.5;
    end
end

for k=1:sense_node    
    for i=1:grid_range_y*span  
        for j=1:grid_range_x*span  
            dist_node_target(i,j,k)=dist(node_x(k),node_y(k),target_x(i,j),target_y(i,j));
        end
    end
end

Eo=0.25;
ETX=50*0.000000001;
ERX=50*0.000000001;
%Transmit Amplifier types
Emp=0.1*0.000000001;
%Data Aggregation Energy
EDA=5*0.000000001;
flag=0;  
do=87;
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

dist_node=zeros(sense_node+1,sense_node+1); 
for i=1:(sense_node+1)
    for j=1:(sense_node+1)
        dist_node(i,j)=dist(S(i).xd,S(i).yd,S(j).xd,S(j).yd);
    end
end

target_covered_for_each_node=zeros(length(target_x(:,1)),length(target_x(1,:)),sense_node); 
for k=1:sense_node 
    for i=1:length(target_x(:,1)) 
        for j=1:length(target_x(1,:))
            if dist_node_target(i,j,k)<=rs
                target_covered_for_each_node(i,j,k)=1;
            end
        end
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

%%%%%%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%% START OF ROUNDS %%%%%%%%%%%%%%%%%%%%%%%%%%
for r=1:1:rmax
  %Operation for epoch
  r
  

    if mod(r, round(1/p) )==0
       for i=1:sense_node  
         S(i).G=0;
         S(i).cl=0;
       end
    end
    for y=1:sense_node
       S(y).sons=0; 
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
%         plot(S(i).xd,S(i).yd,'red .');
        dead=dead+1;
%         hold on;    
    end

    if S(i).E>0
        S(i).type='N';
    end
    S(i).sons=0; 
end
% plot(S(n+1).xd,S(n+1).yd,'x');

coverage_rec(r+1)=eval_coverage(S,target_covered_for_each_node,target_x,sense_node);
    

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
cluster=1; 
alive_node=0;
for i=1:1:sense_node
   if(S(i).E>0)
   temp_rand=rand;    
   alive_node=alive_node+1;
   if ( (S(i).G)<=0) 
       %Election of Cluster Heads    
       if(temp_rand<= (p/(1-p*mod(r,round(1/p)))))
            countCHs=countCHs+1;
            packets_TO_BS=packets_TO_BS+1;            
            S(i).type='C';
            S(i).G=round(1/p)-1;  
            C(cluster).xd=S(i).xd;
            C(cluster).yd=S(i).yd;
            
            distance=dist_node(i,sense_node+1);
            C(cluster).distance=distance;
            C(cluster).id=i;
            X(cluster)=S(i).xd;
            Y(cluster)=S(i).yd;
            cluster=cluster+1;
        end     
    
    end
   end 

end

STATISTICS(r+1).CLUSTERHEADS=cluster-1;
CLUSTERHS(r+1)=cluster-1;
if (alive_node==0 && flag==0 )||(r==rmax && flag==0 ) 
    last_round=r;
    flag=1;
end
%Election of Associated Cluster Head for Normal Nodes
for i=1:1:sense_node
   if ( S(i).type=='N' && S(i).E>0 )
     if(cluster-1>=1) 
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
        S(i).E=S(i).E- (EDA*(packet_bit) + ETX*(packet_bit) + Emp*packet_bit*( min_dis * min_dis)); 
        packets_TO_BS=packets_TO_BS+1; 
     end
 end
end
      PACKETS_TO_CH(r+1)=packets_TO_CH; 
      PACKETS_TO_BS(r+1)=packets_TO_BS;        
          
for h=1:sense_node 
    if S(h).type=='C' && S(h).E>0
       dis=dist_node(h,sense_node+1);
       S(h).E=S(h).E- ( EDA*(packet_bit) + ETX*(packet_bit) + Emp*packet_bit*( dis * dis)); 
    end
end
hold on;

countCHs;
rcountCHs=rcountCHs+countCHs;

for i=1:sense_node
    if S(i).E<0
        S(i).E=0;
    end
end

end

figure(2);
if last_round~=0
    avg_ch=sum(CLUSTERHS)/last_round;
    avg_packets_to_bs=sum(PACKETS_TO_BS)/last_round;
    avg_packets_to_ch=sum(PACKETS_TO_CH)/last_round;
else
    avg_ch=sum(CLUSTERHS)/rmax;
end
plot(sense_node-DEAD(2:(rmax+1)),'b');  
figure(3);
target_num=length(target_x(1,:))*length(target_y(:,1));
plot(coverage_rec(2:(rmax+1))/target_num,'r'); 



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




