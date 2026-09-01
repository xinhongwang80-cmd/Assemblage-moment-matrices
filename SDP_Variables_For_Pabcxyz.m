function [pax,pby,pcz,pabxy,pacxz,pbcyz,pabcxyz] = SDP_Variables_For_Pabcxyz(nx, ny, nz, na, nb, nc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SDP_Variables_For_Pabcxyz.m generates a set of SDP varialbes for
% probability distributions P(a,b,c|x,y,z)
%
% nx, ny, nz : number of measurement settings for Alice, Bob, and Charlie
% na, nb, nc : number of measurement outcomes for Alice, Bob, and Charlie
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pax = cell(na,nx);
pby = cell(nb,ny);
pcz = cell(nc,nz);
pabxy = cell(na,nb,nx,ny);
pacxz = cell(na,nc,nx,nz);
pbcyz = cell(nb,nc,ny,nz);
pabcxyz = cell(na,nb,nc,nx,ny,nz);

for x = 1:nx   
    sum_a_pax{x} = 0;
    for a = 1:na-1
        pax{a,x} = sdpvar(1,1,'hermitian','real');
        sum_a_pax{x} = sum_a_pax{x} + pax{a,x};
    end
    pax{na,x} = 1 - sum_a_pax{x};
end

for y = 1:ny
    sum_b_pby{y} = 0;
    for b = 1:nb-1
        pby{b,y} = sdpvar(1,1,'hermitian','real');
        sum_b_pby{y} = sum_b_pby{y} + pby{b,y};
    end
    pby{nb,y} = 1 - sum_b_pby{y};
end

for z = 1:nz
    sum_c_pcz{z} = 0;
    for c = 1:nc-1
        pcz{c,z} = sdpvar(1,1,'hermitian','real');
        sum_c_pcz{z} = sum_c_pcz{z} + pcz{c,z};
    end
    pcz{nc,z} = 1 - sum_c_pcz{z};
end


for x = 1:nx
    for y = 1:ny
        sum_ab_pabxy{x,y} = 0;
        for a = 1:na-1
            for b = 1:nb-1
                pabxy{a,b,x,y} = sdpvar(1,1,'hermitian','real');
                sum_ab_pabxy{x,y} = sum_ab_pabxy{x,y} + pabxy{a,b,x,y};
            end
        end
    end
end

for x = 1:nx
    for y = 1:ny
        sum_a_pabxy_fix_b{x,y} = 0;
        for a = 1:na-1
            sum_b_pabxy{a,x,y} = 0;
            for b = 1:nb-1
                sum_b_pabxy{a,x,y} = sum_b_pabxy{a,x,y} + pabxy{a,b,x,y};
            end
            pabxy{a,nb,x,y} = pax{a,x} - sum_b_pabxy{a,x,y};
            sum_a_pabxy_fix_b{x,y} = sum_a_pabxy_fix_b{x,y} + pabxy{a,nb,x,y};
        end
    end
end

for x = 1:nx
    for y = 1:ny
        sum_b_pabxy_fix_a{x,y} = 0;
        for b = 1:nb-1
            sum_a_pabxy{b,x,y} = 0;
            for a = 1:na-1
                sum_a_pabxy{b,x,y} = sum_a_pabxy{b,x,y} + pabxy{a,b,x,y};
            end
            pabxy{na,b,x,y} = pby{b,y} - sum_a_pabxy{b,x,y};
            sum_b_pabxy_fix_a{x,y} = sum_b_pabxy_fix_a{x,y} + pabxy{na,b,x,y};
        end
    end
end

for x = 1:nx
    for y = 1:ny
        pabxy{na,nb,x,y} = 1 - sum_ab_pabxy{x,y} - sum_a_pabxy_fix_b{x,y} - sum_b_pabxy_fix_a{x,y};
    end
end


for x = 1:nx
    for z = 1:nz
        sum_ac_pacxz{x,z} = 0;
        for a = 1:na-1
            for c = 1:nc-1
                pacxz{a,c,x,z} = sdpvar(1,1,'hermitian','real');
                sum_ac_pacxz{x,z} = sum_ac_pacxz{x,z} + pacxz{a,c,x,z};
            end
        end
    end
end

for x = 1:nx
    for z = 1:nz
        sum_a_pacxz_fix_c{x,z} = 0;
        for a = 1:na-1
            sum_c_pacxz{a,x,z} = 0;
            for c = 1:nc-1
                sum_c_pacxz{a,x,z} = sum_c_pacxz{a,x,z} + pacxz{a,c,x,z};
            end
            pacxz{a,nc,x,z} = pax{a,x} - sum_c_pacxz{a,x,z};
            sum_a_pacxz_fix_c{x,z} = sum_a_pacxz_fix_c{x,z} + pacxz{a,nc,x,z};
        end
    end
end

for x = 1:nx
    for z = 1:nz
        sum_c_pacxz_fix_a{x,z} = 0;
        for c = 1:nc-1
            sum_a_pacxz{c,x,z} = 0;
            for a = 1:na-1
                sum_a_pacxz{c,x,z} = sum_a_pacxz{c,x,z} + pacxz{a,c,x,z};
            end
            pacxz{na,c,x,z} = pcz{c,z} - sum_a_pacxz{c,x,z};
            sum_c_pacxz_fix_a{x,z} = sum_c_pacxz_fix_a{x,z} + pacxz{na,c,x,z};
        end
    end
end

for x = 1:nx
    for z = 1:nz
        pacxz{na,nc,x,z} = 1 - sum_ac_pacxz{x,z} - sum_a_pacxz_fix_c{x,z} - sum_c_pacxz_fix_a{x,z};
    end
end

for y = 1:ny
    for z = 1:nz
        sum_bc_pbcyz{y,z} = 0;
        for b = 1:nb-1
            for c = 1:nc-1
                pbcyz{b,c,y,z} = sdpvar(1,1,'hermitian','real');
                sum_bc_pbcyz{y,z} = sum_bc_pbcyz{y,z} + pbcyz{b,c,y,z};
            end
        end
    end
end

for y = 1:ny
    for z = 1:nz
        sum_b_pbcyz_fix_c{y,z} = 0;
        for b = 1:nb-1
            sum_c_pbcyz{b,y,z} = 0;
            for c = 1:nc-1
                sum_c_pbcyz{b,y,z} = sum_c_pbcyz{b,y,z} + pbcyz{b,c,y,z};
            end
            pbcyz{b,nc,y,z} = pby{b,y} - sum_c_pbcyz{b,y,z};
            sum_b_pbcyz_fix_c{y,z} = sum_b_pbcyz_fix_c{y,z} + pbcyz{b,nc,y,z};
        end
    end
end

for y = 1:ny
    for z = 1:nz
        sum_c_pbcyz_fix_b{y,z} = 0;
        for c = 1:nc-1
            sum_b_pbcyz{c,y,z} = 0;
            for b = 1:nb-1
                sum_b_pbcyz{c,y,z} = sum_b_pbcyz{c,y,z} + pbcyz{b,c,y,z};
            end
            pbcyz{nb,c,y,z} = pcz{c,z} - sum_b_pbcyz{c,y,z};
            sum_c_pbcyz_fix_b{y,z} = sum_c_pbcyz_fix_b{y,z} + pbcyz{nb,c,y,z};
        end
    end
end

for y = 1:ny
    for z = 1:nz
        pbcyz{nb,nc,y,z} = 1 - sum_bc_pbcyz{y,z} - sum_b_pbcyz_fix_c{y,z} - sum_c_pbcyz_fix_b{y,z};
    end
end

for x = 1:nx
    for y = 1:ny
        for z = 1:nz
            sum_abc_pabcxyz{x,y,z} = 0;
            for a = 1:na-1
                for b = 1:nb-1
                    for c = 1:nc-1
                        pabcxyz{a,b,c,x,y,z} = sdpvar(1,1,'hermitian','real');
                        sum_abc_pabcxyz{x,y,z} = sum_abc_pabcxyz{x,y,z} + pabcxyz{a,b,c,x,y,z};
                    end
                end
            end
        end
    end
end

for x = 1:nx
    for y = 1:ny
        for z = 1:nz
            for b = 1:nb-1
                sum_ac_pabcxyz{b,x,y,z} = 0;
                for c = 1:nc-1
                    sum_a_pabcxyz{b,c,x,y,z} = 0;
                    for a = 1:na-1
                        sum_a_pabcxyz{b,c,x,y,z} = sum_a_pabcxyz{b,c,x,y,z} + pabcxyz{a,b,c,x,y,z};
                        sum_ac_pabcxyz{b,x,y,z} = sum_ac_pabcxyz{b,x,y,z} + pabcxyz{a,b,c,x,y,z};
                    end
                    pabcxyz{na,b,c,x,y,z} = pbcyz{b,c,y,z} - sum_a_pabcxyz{b,c,x,y,z};
                end
                
            end
        end
    end
end

for x = 1:nx
    for y = 1:ny
        for z = 1:nz
            for c = 1:nc-1
                sum_ab_pabcxyz{c,x,y,z} = 0;
                for a = 1:na-1
                    sum_b_pabcxyz{a,c,x,y,z} = 0;
                    for b = 1:nb-1
                        sum_b_pabcxyz{a,c,x,y,z} = sum_b_pabcxyz{a,c,x,y,z} + pabcxyz{a,b,c,x,y,z};
                        sum_ab_pabcxyz{c,x,y,z} = sum_ab_pabcxyz{c,x,y,z} + pabcxyz{a,b,c,x,y,z};
                    end
                    pabcxyz{a,nb,c,x,y,z} = pacxz{a,c,x,z} - sum_b_pabcxyz{a,c,x,y,z};
                end
                
            end
        end
    end
end

for x = 1:nx
    for y = 1:ny
        for z = 1:nz
            for a = 1:na-1
                sum_bc_pabcxyz{a,x,y,z} = 0;
                for b = 1:nb-1
                    sum_c_pabcxyz{a,b,x,y,z} = 0;
                    for c = 1:nc-1
                        sum_c_pabcxyz{a,b,x,y,z} = sum_c_pabcxyz{a,b,x,y,z} + pabcxyz{a,b,c,x,y,z};
                        sum_bc_pabcxyz{a,x,y,z} = sum_bc_pabcxyz{a,x,y,z} + pabcxyz{a,b,c,x,y,z};
                    end
                    pabcxyz{a,b,nc,x,y,z} = pabxy{a,b,x,y} - sum_c_pabcxyz{a,b,x,y,z};
                end
                
            end
        end
    end
end

for x = 1:nx
    for y = 1:ny
        for z = 1:nz
            for a = 1:na-1
                sum_b_pabcxyz_fix_c{a,x,y,z} = 0;
                for b = 1:nb-1
                    sum_b_pabcxyz_fix_c{a,x,y,z} = sum_b_pabcxyz_fix_c{a,x,y,z} + pabcxyz{a,b,nc,x,y,z};
                end
                sum_c_pabcxyz_fix_b{a,x,y,z} = 0;
                for c = 1:nc-1
                    sum_c_pabcxyz_fix_b{a,x,y,z} = sum_c_pabcxyz_fix_b{a,x,y,z} + pabcxyz{a,nb,c,x,y,z};
                end
            end
            
            for b = 1:nb-1
                sum_a_pabcxyz_fix_c{b,x,y,z} = 0;
                for a = 1:na-1
                    sum_a_pabcxyz_fix_c{b,x,y,z} = sum_a_pabcxyz_fix_c{b,x,y,z} + pabcxyz{a,b,nc,x,y,z};
                end
                sum_c_pabcxyz_fix_a{b,x,y,z} = 0;
                for c = 1:nc-1
                    sum_c_pabcxyz_fix_a{b,x,y,z} = sum_c_pabcxyz_fix_a{b,x,y,z} + pabcxyz{na,b,c,x,y,z};
                end
            end
            
            for c = 1:nc-1
                sum_a_pabcxyz_fix_b{c,x,y,z} = 0;
                for a = 1:na-1
                    sum_a_pabcxyz_fix_b{c,x,y,z} = sum_a_pabcxyz_fix_b{c,x,y,z} + pabcxyz{a,nb,c,x,y,z};
                end
                sum_b_pabcxyz_fix_a{c,x,y,z} = 0;
                for b = 1:nb-1
                    sum_b_pabcxyz_fix_a{c,x,y,z} = sum_b_pabcxyz_fix_a{c,x,y,z} + pabcxyz{na,b,c,x,y,z};
                end
            end
            
        end
    end
end



for x = 1:nx
    for y = 1:ny
        for z = 1:nz
            
            for a = 1:na-1
                pabcxyz{a,nb,nc,x,y,z} = pax{a,x} - sum_b_pabcxyz_fix_c{a,x,y,z} - sum_c_pabcxyz_fix_b{a,x,y,z} - sum_bc_pabcxyz{a,x,y,z};
            end
            for b = 1:nb-1
                pabcxyz{na,b,nc,x,y,z} = pby{b,y} - sum_c_pabcxyz_fix_a{b,x,y,z} - sum_a_pabcxyz_fix_c{b,x,y,z} - sum_ac_pabcxyz{b,x,y,z};
            end
            for c = 1:nc-1
                pabcxyz{na,nb,c,x,y,z} = pcz{c,z} - sum_a_pabcxyz_fix_b{c,x,y,z} - sum_b_pabcxyz_fix_a{c,x,y,z} - sum_ab_pabcxyz{c,x,y,z};
            end
            
            sum_bc_pabcxyz_fix_a{x,y,z} = sum_bc_pbcyz{y,z} - sum_abc_pabcxyz{x,y,z};
            sum_ac_pabcxyz_fix_b{x,y,z} = sum_ac_pacxz{x,z} - sum_abc_pabcxyz{x,y,z};
            sum_ab_pabcxyz_fix_c{x,y,z} = sum_ab_pabxy{x,y} - sum_abc_pabcxyz{x,y,z};
            
            sum_a_pabcxyz_fix_bc{x,y,z} = sum_a_pax{x} - sum_ab_pabcxyz_fix_c{x,y,z} - sum_ac_pabcxyz_fix_b{x,y,z} - sum_abc_pabcxyz{x,y,z};
            sum_b_pabcxyz_fix_ac{x,y,z} = sum_b_pby{y} - sum_ab_pabcxyz_fix_c{x,y,z} - sum_bc_pabcxyz_fix_a{x,y,z} - sum_abc_pabcxyz{x,y,z};
            sum_c_pabcxyz_fix_ab{x,y,z} = sum_c_pcz{z} - sum_ac_pabcxyz_fix_b{x,y,z} - sum_bc_pabcxyz_fix_a{x,y,z} - sum_abc_pabcxyz{x,y,z};
            
            pabcxyz{na,nb,nc,x,y,z} = 1 - sum_abc_pabcxyz{x,y,z} - sum_bc_pabcxyz_fix_a{x,y,z} - sum_ac_pabcxyz_fix_b{x,y,z} - sum_ab_pabcxyz_fix_c{x,y,z} + ...
                - sum_a_pabcxyz_fix_bc{x,y,z} - sum_b_pabcxyz_fix_ac{x,y,z} - sum_c_pabcxyz_fix_ab{x,y,z};
        end
    end
end





end