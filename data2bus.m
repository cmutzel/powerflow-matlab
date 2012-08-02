function [YbusRe,data,numBus, lines] = data2bus()

load BusInputData.csv;
%load Test3bus.csv;

data = BusInputData;
%data = Test3bus;
[a,b] = size(data);
%data has format: first n rows list the bus numbers in the problem
% row n+1 holds -999 in first column, marking end input bus #s
% after -999 is a NaN row, next m entries are impedance data 
% n determines the size of Ybus

%Find the first row containing -999
BusIndex = 1;
stop = 0;

%Find the end of the bus numbering
while stop == 0
    if data(BusIndex,1) == -999;
        stop = 1; 
    end
    if (stop == 0);
        BusIndex = BusIndex + 1;
    end
    if BusIndex > a;   % end of file 
        stop = 1;
    end
end

%Read in the bus numbers and create an internal mapping scheme 
minBus = min(data(1:BusIndex-1,1));
maxBus = max(data(1:BusIndex-1,1));
mapS = minBus - 1; % External number = internal number + mapS
numBus = maxBus - mapS;

clear Ybus;
%Initialize a sparse matrix element
Ybus = [maxBus - mapS, maxBus - mapS, 0];


% Find stop point in the data
stop = 0;
ConnEnd = BusIndex+1;
while stop == 0
    if data(ConnEnd, 1) == -999;
        stop = 1;
    end
    if (stop == 0);
        ConnEnd = ConnEnd + 1;
    end
end
lines = zeros(ConnEnd - BusIndex -1 ,2);
%for each connection, update Ybus, also store connections
%lines = zeros(length(BusIndex+1:ConnEnd-1),2);
for n = BusIndex+1:ConnEnd-1
    NodeA = data(n,1) - mapS; NodeB = data(n,2) - mapS;
    lines(n - BusIndex,1) = NodeA;
    lines(n - BusIndex,2) = NodeB;
    % Self terms = diagonal elements
    yab = 1/(data(n,6) + 1i*data(n,7));
    Ybus = sparseAdd(Ybus, NodeA, NodeA, yab +.5i*data(n,8));
    Ybus = sparseAdd(Ybus, NodeB, NodeB, yab +.5i*data(n,8));
    % Off diagonal terms, only fill the upper diagonal
    if NodeA > NodeB
        Ybus = sparseAdd(Ybus, NodeB, NodeA, -yab);
    else
        Ybus = sparseAdd(Ybus, NodeA, NodeB, -yab);
    end
end 

%Construct the real matrix Ybus
YbusRe = zeros(Ybus(end,1:2));

for m=1:length(Ybus)-1;
    YbusRe(Ybus(m,1),Ybus(m,2)) = Ybus(m,3);
    if Ybus(m,1) ~= Ybus(m,2)
        YbusRe(Ybus(m,2),Ybus(m,1)) = Ybus(m,3);
    end
end


end