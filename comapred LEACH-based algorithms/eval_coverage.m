function y=eval_coverage(S,target_covered_for_each_node,target_x,sense_node)
    target_covered=zeros(length(target_x(:,1)),length(target_x(:,1)));
    for i=1:sense_node
        if S(i).E>0 
            target_covered(:,:)=or(target_covered(:,:),target_covered_for_each_node(:,:,i));
        end
    end
    y=length(find(target_covered(:,:)==1));
end