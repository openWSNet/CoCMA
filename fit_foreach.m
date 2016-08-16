function [covered_number,covered_map]=fit_foreach(gene) 
  global sense_node target_x target_covered_for_each_node

  target_covered=zeros(length(target_x(:,1)),length(target_x(:,1)));
    for k=1:sense_node
        if(gene(k)==1)
            target_covered(:,:)=or(target_covered(:,:),target_covered_for_each_node(:,:,k));
        end
    end


    covered_target_count=sum(sum(target_covered(:,:),2));
    
    covered_number=covered_target_count;
    covered_map=target_covered;