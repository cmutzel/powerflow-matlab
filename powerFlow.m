%% POWER FLOW ANALYSIS 

% THIS CODE IS ANALYSES THE COMPLEX POWER FLOW FOR A NETWORK

%Files needed in directory are:
% BusInputData.csv
% createYbus.m -forms Ybus from busInputData
% data2bus.m - stores bus data in new structure for calcs
% findJacob.m - finds jacobian
% sparseAdd.m - used in Ybus
% calcPQ.m - calculations PQ injections using calculated node V and Th

% CODE BY CHRIS MUTZEL
% CHRIS.MUTZEL@GMAIL.COM

% EE197 - POWER SYSTEMS ANALYSIS
% PROFESSOR A. STANKOVIC
% DEPEARTMENT OF ELECTRICAL AND COMPUTER ENGINEERING
% TUFTS UNIVERSITY, MEDFORD, MA

% DECEMEBER 21TH, 2011


clc
clear
format short  % keep output to 4 decimal places
disp('NR Power Flow')
disp('Code by Chris Mutzel')
disp('chris.mutze@gmail.com')
Ttotal = tic;

%load Ybus
[Ybus, data, N, lines] = data2bus(); 
save Ybus;

%Create new storage for bustype, P, Q, V, theta
busData = zeros(N,5);           
for n=1:N    %populate dataBus
    busData(data(n,1),1) = data(n,3);                   % bus type
    busData(data(n,1),2) = (data(n,8) - data(n,6))/100; % Pspec [MW]
    busData(data(n,1),3) = (data(n,9) - data(n,7))/100; % Qspec [MVAR]
    busData(data(n,1),4) = data(n,4); % V 
end                                 
disp(' '); disp('Specified Injections'); disp('    #         Pspec     Qspec');
disp([(1:N)' busData(:,2) busData(:,3)])
busData(1,4) = 1.0; %Slack voltage is one, angle zero

numPV=0;
for i = 1:N  %find number of each bus type
   if (busData(i,1) == 2)
      numPV = numPV + 1; %m
   end
end
m = numPV; % N-m-1 are PQ buses

%Intialize N-R Variables
ep  = .00000001;                %stopping criteria for bus power mismatches
STOP = 0;                   %set to one when all deltaPQ < ep
 
X = [zeros(N-1,1); 1*ones(N-m-1,1)];    %flat start
deltaPQ = zeros(2*N-m-2,1);             %Mismatch vector
j = 0;                                  %iteration count

Traph = tic;                       %Start Timer
while STOP == 0;
    if j == 0;
        disp('######################################################')
        disp('Begin Newton Raphson')
        disp('Flat Start->')
        X
    end
    
  
    itNum = int2str(j + 1); 
    disp('######################################################')
    disp(' '); disp(strcat('Begin Iteration: ',itNum)); disp(' ')
    
    % Update V in dataBus for all PQ , use voltages from here for PQcalc
    c = 1;
    for n=1:N
        if busData(n,1) == 0; % if the bus is a PQ bus
           busData(n,4) = X(N+c-1); % store its V along with knowns for calc
        c = c+1;
        end
    end
    
    %Update thetaBus for all PQ, use for PQcalc
     c = 1;
     for n=1:N
         if n == 1;
            busData(n,5) = 0;
         end
         if (busData(n,1) == 2 || busData(n,1) == 0);
            busData(n,5) = X(c); 
         c = c+1;
         end
         
     end
         
    % ###################################################### 
    %Calculate PQcalc based on V and theta
    [Pcalc, Qcalc] = calcPQ(Ybus, busData, N);
    %disp('    #      Pspec    Qspec')
    %disp([(1:N)' busData(:,2) busData(:,3)])
    %disp('    #      Pcalc     Qcalc')
    %disp([(1:N)' Pcalc Qcalc])
   
    % ######################################################
    %Calculate bus power mismatches
    
    %P
    deltaPQ(1:N-1,1) = busData(2:N, 2) - Pcalc(2:N,1);
    
    %Q
    bc = 0; %counter for index in delta PQ
    for n=1:N
        if (busData(n,1) == 0) % only for PQ buses
            bc = bc + 1;
            deltaPQ(N+bc-1,1) = busData(n,3) - Qcalc(n,1); 
        end
    end
    %clear bc
    %disp('Mismatch vector =')
    %disp(deltaPQ)
    
    % ######################################################
    % Check for convergence
    error = norm(deltaPQ,2);
    disp('Error, given by 2-norm of Mismatch vector')
    disp(error)
          
    if (error < ep);
        STOP = 1;
        disp('Power Flow Calculation stopped')
        disp('Error Threshold Reached')
        disp(' ')
    end
    
    % ######################################################
    % New iteration if error over threshold
    if STOP==0;
        
        %Find new Jacobian
        %if j <= 2;
            Jacob = findJacob(Ybus, Pcalc, Qcalc, busData, N, m);
        %end
            
    % Solve the mismatch equations
    % Invert the Jacobian
        invJ = inv(Jacob);
    % Solve for deltaX 
        dX = invJ*deltaPQ;
        
    %Update Theta and V    
        X = X + dX;
    end
    j = j+1;
    %if (STOP == 0)
    %    disp('    #         V       Theta')
    %    disp([(1:N)' busData(:,4) busData(:,5)*360/(2*pi)])
    %    disp(strcat('Finished Iteration:', itNum))
    %end
    disp(' ')
    if j >= 10
        STOP = 1; disp(' ');
        disp('Power Flow Calculation stopped')
        disp('Maximum Iterations (10) Reached'); disp(' ')
    end
      
end   %NR iteration
Traph = toc(Traph); %Stop timer

% ######################################################
% Newton Raphson Finished
% Calculate line flows
Zbus = inv(Ybus);
lineFlows = [lines zeros(length(lines), 1)];
for n = 1:length(lineFlows)
    lineFlows(n,3) =  (exp(1i*angle(Zbus(lineFlows(n,1),lineFlows(n,2))))/abs(Zbus(lineFlows(n,1),lineFlows(n,2))))*...
        (abs(busData(lineFlows(n,1),4))^2 - abs(busData(lineFlows(n,1),4))*abs(busData(lineFlows(n,2),4))*...
         exp(1i*(busData(lineFlows(n,1),5) - busData(lineFlows(n,2),5))));      
end



Ttotal = toc(Ttotal);
disp(['Final Solution after ', int2str(j), ' iterations is:'])
disp('    #         V       Theta')
disp([(1:N)' busData(:,4) busData(:,5)*360/(2*pi)])
disp('Line Flows are:')
disp('  Bus n      to      Bus m              Snm (per unit)')
disp(lineFlows) 
disp('Computation Time for NR-iteration was:')
disp([num2str(Traph), ' seconds'])
disp('Computation Time for program was:')
disp([num2str(Ttotal), ' seconds'])



