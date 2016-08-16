function fitness=GA_LinRank(object, select_pressure)
%% function fitness=GA_LinRank(object, select_pressure)
%% Change objective function values to fitness function ones

% PenChen Chou.  2002-10-25
% Initiallzations
Nind=length(object);
% SP is at most to be 2.
SP=max([0, min([2,select_pressure])]);

% Sort into right position, take the index.
[temp, inx]=sort(object);

% Linear ranking with select_pressure.
% The last position gets the best fitness
for pos=1:Nind
    fitness(pos)=2-SP+2*(SP-1)*(pos-1)/(Nind-1);
end

% Set the fitness back to the right position.
fitness=[fitness(inx)]';
%------------  END ------------------------------------