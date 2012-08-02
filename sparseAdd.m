function newSparseM = sparseAdd(M, row, col, v)
%Function for manipulating sparse matrices
% Input:
%       M is a sparse matrix object 
%       row and col are location in a real matrix M to which v should be added
% Output: 
%       new sparse matrix element O

% M is ordered by row then by column, third column holds value for row and
% col.  Last row has holds dimensions of the matrix object

    
% See if the element exists
loc = 0;
found = 0;       % Loc is filled with new position for an element if found  = 0
c1_index = 1;   % Starting row to check M for existing element
endloop = 0;    % 
[a,b] = size(M); % store dimensions of M

while (endloop ~= 1);    
    if M(c1_index, 1) == row;
        if M(c1_index, 2) == col;   % the element exists, stores its location
            loc = c1_index;
            found = 1;
            endloop = 1;
        end
        if M(c1_index, 2) > col;    % element for same row but higher column   
            loc = c1_index;
            found = 0;            
            endloop = 1;
        end
    end
    
    if M(c1_index, 1) > row;
            loc = c1_index;
            found = 0;
            endloop = 1; 
    end
    c1_index = c1_index + 1;
    if (c1_index > a)  %Check that we are not out of bounds of the array
       found = 0;
       endloop = 1;
    end
end

% if it does not exist, create it
% array has same form up to newly added element, then new element, then
% remaining contents
if found == 0;
    newSparseM = [M(1:loc-1,:); [row, col, v]; M(loc:end,:)];
end

%If does exist, v is added to element value
if found == 1;
    M(loc,3) = M(loc,3) + v;
    newSparseM = M;
end

%insert code here to delete row if value is now zero, not important for
%Ybus since all values remain non-zero once created.
end




