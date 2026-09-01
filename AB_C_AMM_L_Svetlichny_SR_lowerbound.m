clear all;
tic;

bellvalue = 4*sqrt(2);

nx = 2;
ny = 2;
nz = 2;
na = 2;
nb = 2;
nc = 2;

Bell = linspace(4,bellvalue,20);
minSR = zeros(4,20);

nalambda = na^nx;
nblambda = nb^ny;

for level = 1:4
for num = 1:20


S = AMM_proj_gen_xlevel_seq_C(nz,nc,level);

for i = 1:length(S)
    
    Sdag(i) = AMM_proj_adjoint_poly_ABC(S(i));

end

gamma_str_complex = AMM_proj_moment_matrix_string_C(S,Sdag,'complex');
gamma_str_real = AMM_proj_moment_matrix_string_C(S,Sdag,'real');

uni_mono = unique(gamma_str_real);


for a = 1:na
    for b = 1:nb
        for x = 1:nx
            for y = 1:ny

                gamma_SDP{a,b,x,y} = zeros(length(gamma_str_complex));
                gamma_SDP{a,b,x,y} = gamma_SDP{a,b,x,y} + (gamma_str_complex == string('1'));

            end
        end
    end
end



uni_mono(uni_mono == string('1')) = [];
uni_mono(uni_mono == string('0')) = [];


for i = 1:nalambda
    for j = 1:nblambda

        gamma_SDP_lambda{i,j} = zeros(length(gamma_str_complex));
    end
end

for idx = 1:length(uni_mono)

    tfMatrix = (gamma_str_real == uni_mono(idx));
    
    for a = 1:na
        for b = 1:nb
            for x = 1:nx
                for y = 1:ny
                    u{idx,a,b,x,y} = sdpvar(1,1,'hermitian','real');
                    gamma_SDP{a,b,x,y} = gamma_SDP{a,b,x,y} + u{idx,a,b,x,y}*tfMatrix;
                end
            end
        end
    end
end


for i = 1:nalambda
    for j = 1:nblambda
        for idx = 1:length(uni_mono)
            tfMatrix = (gamma_str_real == uni_mono(idx));
            v{idx,i,j} = sdpvar(1,1,'hermitian','real');
            gamma_SDP_lambda{i,j} = gamma_SDP_lambda{i,j} + v{idx,i,j}*tfMatrix;
        end
        v{idx+1,i,j} = sdpvar(1,1,'hermitian','real');
        gamma_SDP_lambda{i,j}(1,1) = v{idx+1,i,j};
    end
end

sum_gamma_SDP_lambda_11 = 0;

for i = 1:nalambda
    for j = 1:nblambda

        sum_gamma_SDP_lambda_11 = sum_gamma_SDP_lambda_11 + gamma_SDP_lambda{i,j}(1,1);

    end
end

DA = detPax(nx,na);
DB = detPby(ny,nb);

for a = 1:na
    for b = 1:nb
        for x = 1:nx
            for y = 1:ny

                sum_D_gamma_sigma_lambda{a,b,x,y} = 0;

                for i = 1:nalambda
                    for j = 1:nblambda
                        
                        sum_D_gamma_sigma_lambda{a,b,x,y} = sum_D_gamma_sigma_lambda{a,b,x,y} + ...
                            DA(i,a,x).*DB(j,b,y).*gamma_SDP_lambda{i,j};

                    end
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


SE = E{1,1,1} + E{1,1,2} + E{1,2,1} - E{1,2,2} + E{2,1,1} - E{2,1,2} - E{2,2,1} - E{2,2,2};


constr = [];


for z = 1:nz
    for c = 1:nc-1
        str_Ccz = string(strcat('C_',num2str(c),'|',num2str(z)));

        for a = 1:na
            for b = 1:nb
                for x = 1:nx
                    for y = 1:ny

                        constr = [constr, pabcxyz{a,b,c,x,y,z} == u{uni_mono==str_Ccz,a,b,x,y}];
                    end
                end
            end
        end
    end
end

for x = 1:nx
    for y = 1:ny
        sum_a_gamma_SDP{x,y} = 0;

        for a = 1:na
            for b = 1:nb
                gamma_SDP{a,b,x,y}(1,1) = pabxy{a,b,x,y};
                sum_a_gamma_SDP{x,y} = sum_a_gamma_SDP{x,y} + gamma_SDP{a,b,x,y};
            end
        end
    end
end

for a = 1:na
    for b = 1:nb
        for x = 1:nx
            for y = 1:ny
                gamma_SDP{a,b,x,y}(1,1) = pabxy{a,b,x,y};
                constr = [constr, gamma_SDP{a,b,x,y} >= 0];
                constr = [constr, sum_D_gamma_sigma_lambda{a,b,x,y} >= gamma_SDP{a,b,x,y}];
            end
        end
    end
end


for x = 1:nx
    for y = 1:ny
        for b = 1:nb
            sum_gamma_SDP_a{b,x,y}=0;

            for a = 1:na

                sum_gamma_SDP_a{b,x,y} = sum_gamma_SDP_a{b,x,y} + gamma_SDP{a,b,x,y};
            end
        end
    end
end


for x = 1:nx
    for y = 1:ny
        for a = 1:na
            sum_gamma_SDP_b{a,x,y}=0;

            for b = 1:nb

                sum_gamma_SDP_b{a,x,y} = sum_gamma_SDP_b{a,x,y} + gamma_SDP{a,b,x,y};
        
            end
        end
    end
end


for i = 1:nx
    for b = 1:nb
        for x = 1:nx
            for y = 1:ny
                
                if i ~= x 
                    constr = [constr, sum_gamma_SDP_a{b,i,y} == sum_gamma_SDP_a{b,x,y}];
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
                    constr = [constr, sum_gamma_SDP_b{a,x,j} == sum_gamma_SDP_b{a,x,y}];
                end
            end
        end
    end
end



for i = 1:nalambda
    for j = 1:nblambda
        constr = [constr, gamma_SDP_lambda{i,j} >= 0];
    end
end

constr = [constr, SE == Bell(num)];

sol = solvesdp(constr, sum_gamma_SDP_lambda_11);
minSR(level,num) = double(sum_gamma_SDP_lambda_11)-1;

sol
toc

end
end


figure;
ax=gca;
set(gcf, 'Color', 'w');
set(ax,'FontSize',12);

axis([4 4*sqrt(2) 0 0.3333]);
box on;
hold on;

plot(Bell,minSR(1,:),'b-','LineWidth',1.5);
plot(Bell,minSR(2,:),'r--','LineWidth',1.5);
plot(Bell,minSR(3,:),'k-.','LineWidth',1.5);

xlabel('Svetlichny inequality','FontSize',15,'Interpreter','latex');
ylabel('minimal steering robustness','FontSize',15,'Interpreter','latex');

legend({'$\ell=1$','$\ell=2$','$\ell=3$'},'FontSize',15,'Interpreter','latex','Location','northwest');

