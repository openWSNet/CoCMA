figure(11);
x=[1 20 50 70 80 90 100];
% the number of dead node : [1 80 200 280 320 360 400] 
y1=[100 750 1950 2600 2900 3400 4050];
y2=[96 859 2075 2912 3395 4090 5066]; % determine by yourself
y3=[96 250 300 350 400 500 550];
y4=[250 1050 1120 1120 1140 1140 1140 ];
plot(x,y2,x,y1,x,y4,x,y3);
xlabel('Percentage of node death (%)','fontsize',14);
ylabel('Number of rounds','fontsize',14);