function seq = AMM_proj_gen_xlevel_seq_C(nz,nc,level)
%%%%

% AMM_proj_gen_xlevel_seq.m generates the sequence of Charlie's projector
% nx, ny, nz : Number of measurement settings for Alice, Bob, and Charlie.
% na, nb, nc : Number of outcomes per measurement for each party.
% level: moment relaxation level

%%%%
    seq_S_C = [];
    for z = 1:nz
        for c = 1:nc-1

            seq_S_C = [seq_S_C string(strcat('C_',num2str(c),'|',num2str(z)))];
        end
    end

    if level == 1
        
        seqC = [string('Id') seq_S_C];
    else

        seqC = [string('Id') seq_S_C];
        SCcal{1} = seq_S_C;
        
        S = seq_S_C;
        
        for k = 1:level-1
            SCcal{k+1} = [];
            for i = 1:length(SCcal{k})
                for j = 1:length(S)
                    Sdet = split(SCcal{k}(i),'*');
                    if Sdet{length(Sdet)} ~= S(j)

                        Snew = SCcal{k}(i) + '*' + S(j);
                        SCcal{k+1} = [SCcal{k+1}, Snew];
                        seqC = [seqC Snew];
                    end
                end
            end
            
            
        end
    end

    seq = seqC;

end