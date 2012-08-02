function [Pcalc,Qcalc] = calcPQ(Ybus, busData, N)
% Summary: Calculates P and Q injections at each node
% 
% Inputs:  
%   scalar N = numNonSlack is # of non slack busses in the power flow
%   X is column vector with 1st N - 1  entries for non Slack buses [theta]
%                           then N - m -1  entries for PQ buses [voltage]
%
% Outputs:  
%   PQ injections at all buses

%Create storage for output
Pcalc = zeros(N,1); %top row will refer to slack bus
Qcalc = zeros(N,1);

%P
for i=1:N
    for k=1:N
        Pcalc(i,1) = Pcalc(i,1) + busData(i,4)*busData(k,4)*(...
            real(Ybus(i,k))*cos(busData(i,5)-busData(k,5)) +...
            imag(Ybus(i,k))*sin(busData(i,5)-busData(k,5)));
    end
end


for i=1:N
    for k=1:N
        Qcalc(i,1) = Qcalc(i,1) + busData(i,4)*busData(k,4)*(...
            real(Ybus(i,k))*sin(busData(i,5)-busData(k,5)) -...
            imag(Ybus(i,k))*cos(busData(i,5)-busData(k,5)));
    end
end

end

