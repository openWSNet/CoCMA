[x,y]=meshgrid(0:0.01:1,0:0.01:1);
z=(x.^2).*1-(y.^0.5).*1;
mesh(x,y,z);