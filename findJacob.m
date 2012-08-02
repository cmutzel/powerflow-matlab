function [newJacob] = findJacob2(Ybus,Pcalc, Qcalc, busData, N, m)
% Summary: Finds the Jacobian for the power flow problem given
%
% Inputs:
%   scalar N = numNonSlack is # of non slack busses in the power flow
%   X is column vector with 1st N entries theta for non Slack buses
%                           then N entries V for non Slack buses
%
% Outputs:
%   Jacob 

%Create storage for output
%square matrix with n = N-1 + N - m - 1
Jacob = zeros(2*N,2*N);


sizeJ = 2*N;

for Jr=1:2
    for Jc=1:2
        
        %J11
        if (Jr == 1 && Jc == 1)
            for SJr=1:N
                
                for SJc=1:N
                   
                    if (SJr == SJc) % diagonal elements 
                        Jacob(SJr,SJr) = -Qcalc(SJc,1) - imag(Ybus(SJc,SJc))*abs(busData(SJc,4))^2;
                                
                    else            % off-diagonal elements
                        Jacob(SJr, SJc) = abs(busData(SJr,4))*abs(busData(SJc,4))*(...
                            real(Ybus(SJr,SJc))*sin(busData(SJr,5)-busData(SJc,5))-...
                            imag(Ybus(SJr,SJc))*cos(busData(SJr,5)-busData(SJc,5)));
                        
                    end
                end
            end
            
        end
        
        %J12
        if (Jr == 1 && Jc == 2)
            for SJr=1:N
                for SJc=1:N
                    if SJr == SJc % diagonal elements
                        Jacob(SJr, N + SJc) = Pcalc(SJc)/abs(busData(SJc,4)) +...
                            real(Ybus(SJc,SJc))*abs(busData(SJc,4));
                    else          % off-diagaonal elements
                        Jacob(SJr, N + SJc) = abs(busData(SJr,4))*(...
                            real(Ybus(SJr, SJc))*cos(busData(SJr,5)-busData(SJc, 5))+...
                            imag(Ybus(SJr, SJc))*sin(busData(SJr,5)-busData(SJc, 5)));
                    end
                end
            end
        end
        
        
        %J21
        if (Jr == 2 && Jc == 1)
            for SJr=1:N
                for SJc=1:N
                    if SJr == SJc % diagonal elements
                        Jacob(N + SJr, SJc) = Pcalc(SJc) - real(Ybus(SJc, SJc))*abs(busData(SJc, 4))^2;
                    else          % off-diagonal elements
                        Jacob(N + SJr, SJc) = -abs(busData(SJr,4))*abs(busData(SJc,4))*(...
                            real(Ybus(SJr, SJc))*cos(busData(SJr,5) - busData(SJc,5)) +... 
                            imag(Ybus(SJr, SJc))*sin(busData(SJr,5) - busData(SJc,5)));
                    end
                end
            end
        end
        
        %J22
        if (Jr == 2 && Jc == 2)
            for SJr=1:N
                for SJc=1:N
                    if SJr == SJc % diagonal elements
                        Jacob(N + SJr, N +  SJc) = Qcalc(SJc,1)/abs(busData(SJc,4)) - ...
                            imag(Ybus(SJc, SJc))*abs(busData(SJc,4));
                    else          %off-diagonal elements
                        Jacob(N + SJr, N +  SJc) = abs(busData(SJr, 4))*(...
                            real(Ybus(SJr, SJc))*sin(busData(SJr,5)-busData(SJc,5))-...
                            imag(Ybus(SJr, SJc))*cos(busData(SJr,5)-busData(SJc,5)));
                            
                    end
                end
            end
        end
        
    end
end

% ######################################################
%Just found the Jacob for all buses, all only for PQ + PV, PV
newJacob = zeros(2*N-m-2,2*N-m-2);

%J11
for c = 2:N
    newJacob(1:N-1,c-1) = Jacob(2:N, c);
end

%J12
a=0;
for b = 1:N
    if (busData(b,1) == 0)
        a = a + 1;
        newJacob(1:N-1, N-1+a) = Jacob(2:N,b+N);
    end
end

%J21
a=0;
for b=1:N
    if (busData(b,1) == 0)
        a = a + 1;
        newJacob(N-1+a, 1:N-1) = Jacob(b+N, 2:N);
    end
end

%J22
%nest for loops do same thing as above

a=0;
for b=1:N
    if (busData(b,1) == 0)
        a=a+1;
        c=0;
        for d=1:N
            if (busData(d,1) == 0)
                c=c+1;
                newJacob(N-1+a, N-1+c) = Jacob(b+N,d+N);
            end
        end
    end
end
end

