function z=energy_analysis(x_range,y_range,node,node_size)
    x=zeros(y_range,x_range);
    y=zeros(y_range,x_range);
    z=zeros(y_range,x_range,node_size);
    z_sum=zeros(y_range,x_range);
    for p=1:y_range
        for q=1:x_range
           x(p,q)=q;
           y(p,q)=p;
        end
    end
    
    for k=1:node_size
        for i=1:1:y_range
            for j=1:x_range
                z(i,j,k)=node(k).E.*exp(-0.03.*((j-node(k).xd).^2+(i-node(k).yd).^2));
            end
        end
    end
    
    for g=1:20
        margin_top_node(g).x=g*5-5;
        margin_top_node(g).y=100;
        margin_top_node(g).E=node(380+g).E;
        margin_right_node(g).x=100;
        margin_right_node(g).y=g*5-5;
        margin_right_node(g).E=node(g*20).E;
    end
    for g=1:20
        for i=1:1:y_range
            for j=1:x_range
                z_margin_top(i,j,g)=margin_top_node(g).E.*exp(-0.03.*((j-margin_top_node(g).x).^2+(i-margin_top_node(g).y).^2));
                z_margin_right(i,j,g)=margin_right_node(g).E.*exp(-0.03.*((j-margin_right_node(g).x).^2+(i-margin_right_node(g).y).^2));
            end
        end     
    end
    
    for m=1:node_size
        z_sum=z_sum+z(:,:,m);
    end
    
    for g=1:20
        z_sum=z_sum+z_margin_top(:,:,g)+z_margin_right(:,:,g);
    end
    
    figure(10);
    axis equal;
    surf(x,y,z_sum);
    caxis([0 1])
    colorbar;
    shading interp;
    view(0,90);
    colormap('jet');
end