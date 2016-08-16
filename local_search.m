function y=local_search(sensor_select) % do it for each generation
   pop=sensor_select; % get population
        for i=1:length(pop(:,1)) % for each chromosome
            active_node=find(pop(i,:)==1); 
            for j=1:length(active_node)  % find the active nodes and try to deactivate them §ä¥Xactive_node, then reevaluate the new chromosome
                [f_pre,temp1]=fit_foreach(pop(i,:)); % evaluate the chromosome
                pop(i,active_node(j))=0;
                [f_fit,temp2]=fit_foreach(pop(i,:));
                if f_pre>f_fit
                    pop(i,active_node(j))=1;
                end         
            end
        end          
y=pop;        