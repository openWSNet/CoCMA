function gene_output=local_search_for_alive_neighbor(dead_node,gene_input,node,dist_node)
% dead_node : the index of dead node
% gene_input : chromosome input (including the dead node = 0)
% neighbor : record of neighbors
% neighobr_gene : record of which neighbors to be awakened (=1)
% error_range : threshold to 喚醒鄰居的參考值 預設1 (假使不喚醒該鄰居所造成的覆蓋損失大於error_range，則喚醒該鄰居)(跟全部喚醒的覆蓋情況做比較)
global  sense_node sense_range target_covered_for_each_node 
trans_range=17.675; % transmission distance
neighbor=[];
for i=1:sense_node
    if i~=dead_node      
        if dist_node(i,dead_node)<trans_range && gene_input(i)==0 && node(i).E>0 && node(i).type~='D'
            neighbor(length(neighbor)+1)=i;
        end
    end
end

[temp1,no_dead_node_covered_map]=fit_foreach(gene_input);
gene_tmp=gene_input;
gene_tmp(dead_node)=1;
[temp2,include_dead_node_covered_map]=fit_foreach(gene_tmp);
remaining_no_covered_target=xor(no_dead_node_covered_map,include_dead_node_covered_map); % bare targets caused by the dead nodes

% optimal_gene=zeros(1,2);
lost_target_num=length(find(remaining_no_covered_target==1));
optimal_covered_num=0;
optimal_rec=[];
fprintf('\nlossing no.%d',dead_node);
fprintf('\n # of neighbors=%d, # of lost targets=%d',length(neighbor),lost_target_num);

if isempty(neighbor)~=1 
    candidate_node=[];
    for k=1:length(neighbor)
       tmp=and(target_covered_for_each_node(:,:,neighbor(k)),remaining_no_covered_target);
       sum_tmp=sum(sum(tmp,1));
       if  sum_tmp>=1
           candidate_node=[candidate_node neighbor(k)];
       end
    end
%     fprintf('\nthe neighbors covering the lost targets\n');
%     disp(candidate_node);
    
    if length(candidate_node)>4  % limit the max. number of neighbor to be awakened
        cbn_num=4;
    else
        cbn_num=length(candidate_node);
    end
    
    for i=1:cbn_num
        disp(i);
        choice=[];     
        choice=combntns(candidate_node,i);  %% any combinations       
        for j=1:length(choice(:,1))
            count_num=evaluate_with_remaining_area(choice(j,:),remaining_no_covered_target);
%             disp('combination');     
%             disp(choice(j,:));
%             disp('count_num=');     
%             disp(count_num);
            if count_num>optimal_covered_num
                optimal_covered_num=count_num;
                optimal_rec=choice(j,:);
            end
        end
    end
end

gene_output=gene_input;
fprintf('\nwake=');
if optimal_covered_num~=0 
    for m=1:length(optimal_rec)
        gene_output(optimal_rec(m))=1; 
        fprintf('%d ',optimal_rec(m));
    end
end
end


