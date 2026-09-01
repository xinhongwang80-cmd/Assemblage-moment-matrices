clear all;
tic;

na = 2;
nb = 2;
nc = 2;
nx = 2;
ny = 2;
nz = 2;
level = 2;


S = AMM_proj_gen_xlevel_seq_C(nz,nc,level);

for i = 1:length(S)
    
    Sdag(i) = AMM_proj_adjoint_poly_ABC(S(i));

end

gamma_str = AMM_proj_moment_matrix_string_C(S,Sdag,'real');


uni_mono = unique(gamma_str);

for x = 1:nx
    for y = 1:ny
        for a = 1:na
            for b = 1:nb
                gamma_SDP{a,b,x,y} = zeros(length(gamma_str),length(gamma_str));
            end
        end
    end
end


uni_mono(uni_mono==string('0'))=[];

for idx = 1:length(uni_mono)
    
    tfMatrix = (gamma_str == uni_mono(idx)); 

    for a = 1:na
        for b = 1:nb
            for x = 1:nx
                for y = 1:ny
                    u{idx,a,b,x,y} = sdpvar(1,1,'hermitian','real');
                    gamma_SDP{a,b,x,y} = gamma_SDP{a,b,x,y} + u{idx,a,b,x,y}.*tfMatrix;
                end
            end
        end
    end

end

[pax,pby,pcz,pabxy,pacxz,pbcyz,pabcxyz] = SDP_Variables_For_Pabcxyz(nx, ny, nz, na, nb, nc);


for x = 1:nx
    for y = 1:ny
        for z = 1:nz
            E{x,y,z} = pabcxyz{1,1,1,x,y,z} + pabcxyz{1,2,2,x,y,z} + pabcxyz{2,1,2,x,y,z} + ...
                       pabcxyz{2,2,1,x,y,z} - pabcxyz{1,1,2,x,y,z} - pabcxyz{1,2,1,x,y,z} - ...
                       pabcxyz{2,1,1,x,y,z} - pabcxyz{2,2,2,x,y,z};
        
        end
    end
end


SE = E{1,1,1} - E{1,2,2} - E{2,1,2} - E{2,2,1};

constr = [];


for z = 1:nz
    for c = 1:(nc-1)
        str_Ccz = string(strcat('C_',num2str(c),'|',num2str(z)));
        
        for a = 1:na
            for x = 1:nx
                for b = 1:nb
                    for y = 1:ny
                        constr = [constr, pabcxyz{a,b,c,x,y,z} == u{uni_mono==str_Ccz,a,b,x,y}];

                    end
                end
            end
        end
        
    end
end

for a = 1:na
    for b = 1:nb
        for x = 1:nx
            for y = 1:ny

                constr = [constr, gamma_SDP{a,b,x,y}(1,1) == pabxy{a,b,x,y}];
                constr = [constr, pabxy{a,b,x,y} >= 0];
            end
        end
    end
end


for x = 1:nx
    for y = 1:ny
        for b = 1:nb
            sum_gamma_a{b,x,y} = 0;
            for a = 1:na

                sum_gamma_a{b,x,y} = sum_gamma_a{b,x,y} + gamma_SDP{a,b,x,y};
                
            end
        end
    end
end


for x = 1:nx
    for y = 1:ny
        for a = 1:na
            sum_gamma_b{a,x,y} = 0;
            for b = 1:nb
              
                sum_gamma_b{a,x,y} = sum_gamma_b{a,x,y} + gamma_SDP{a,b,x,y};
                
            end
        end
    end
end


for x = 1:nx
    for y = 1:ny
        for a = 1:na
            for b = 1:nb

                constr = [constr, gamma_SDP{a,b,x,y} >= 0];                
            
            end
        end
    end
end


for i = 1:nx
    for b = 1:nb
        for x = 1:nx
            for y = 1:ny
                
                if i ~= x 
                    constr = [constr, sum_gamma_a{b,i,y} == sum_gamma_a{b,x,y}];
                end
            end
        end
    end
end

for j = 1:ny
    for a = 1:na
        for x = 1:nx
            for y = 1:ny
                
                if j ~= y 
                    constr = [constr, sum_gamma_b{a,x,j} == sum_gamma_b{a,x,y}];
                end
            end
        end
    end
end


sol = solvesdp(constr , -SE);
sol
value_SE = double(SE)


toc