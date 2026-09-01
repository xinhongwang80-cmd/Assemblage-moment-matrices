function anti_comm_result_full = AMM_proj_adjoint_poly_ABC(string_of_poly)
%%%%%
% AMM_proj_adjoint_poly_ABC.m Compute the adjoint of a tripartite projector polynomial.
% e.g. transform "A_1|1*A_1|2*B_1|1" into "A_1|2*A_1|1*B_1|1"
% string_of_poly should be in the string form
%%%%%

poly_split = strsplit(string_of_poly,{'+','-'});

char_of_poly = char(string_of_poly);
plus_and_minus = char_of_poly(or(char_of_poly=='+',char_of_poly=='-'));

if poly_split(1)==string()
    poly_split(1) = [];
end

if length(poly_split) == length(plus_and_minus)+1
    % this means the leading coefficient is +, such as
    % 1+A1-A2 or A2-A1B2+B1
    
    plus_and_minus = ['+',plus_and_minus];
    
end

anti_comm_result_full = string();
for i = 1:length(poly_split)
    
gamma = poly_split(i);

if gamma == string('1')
    
    anti_comm_result(i) = string('1');
    
elseif gamma == string('-1')
    
    anti_comm_result(i) = string('-1');
    
else
    
    gamma_sep = strsplit(gamma,'*');
    
    if length(gamma_sep)==1
        anti_comm_result(i) = gamma_sep;
    else
        gamma_A = gamma_sep(~cellfun('isempty', strfind(gamma_sep,'A')));
        % this is to only keep projectors of A
        gamma_B = gamma_sep(~cellfun('isempty', strfind(gamma_sep,'B')));
        % this is to only keep projectors of B
        gamma_C = gamma_sep(~cellfun('isempty', strfind(gamma_sep,'C')));
        if ~isempty(gamma_A) % length of gamma_A =\= 0
            gamma_A_dag = flip(gamma_A);
            gamma_A_dag = strjoin(gamma_A_dag,'*');
        elseif isempty(gamma_A)
            gamma_A_dag = [];
        end
        
        if ~isempty(gamma_B) % length of gamma_B =\= 0
            gamma_B_dag = flip(gamma_B);
            gamma_B_dag = strjoin(gamma_B_dag,'*');
        elseif isempty(gamma_B)
            gamma_B_dag = [];
        end

        if ~isempty(gamma_C) % length of gamma_A =\= 0
            gamma_C_dag = flip(gamma_C);
            gamma_C_dag = strjoin(gamma_C_dag,'*');
        elseif isempty(gamma_C)
            gamma_C_dag = [];
        end
        
        anti_comm_result(i) = strjoin(string([gamma_A_dag gamma_B_dag gamma_C_dag]),'*');
    end
    
end

anti_comm_result_full = strjoin([anti_comm_result_full...
        strjoin([string(plus_and_minus(i)) anti_comm_result(i)],'')],'');

end

anti_comm_result_full_char = char(anti_comm_result_full);
if anti_comm_result_full_char(1) == '+'
    
    anti_comm_result_full_char(1) = [];
    anti_comm_result_full = string(anti_comm_result_full_char);
    
end


end