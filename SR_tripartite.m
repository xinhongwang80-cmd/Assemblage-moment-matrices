function [y1,y2] = SR_tripartite(assemblage)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SR_AB.m Compute the steering robustness of the given assemblage
% with respect to the local LHS model, where Alice and Bob steer Charlie.
% assemblage: 
%   the assemblage created by a set of Alice and Bob's POVMs 
%   shape: [d, d, na, nb, nx, ny]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[d, ~, na, nb, nx, ny] = size(assemblage);

na_lambda = na^nx;
nb_lambda = nb^ny;
nlambda   = na_lambda * nb_lambda;

% hidden states sigma_lambda
for lambda = 1:nlambda
    s{lambda} = sdpvar(d,d,'hermitian','complex'); % σ_λ
end

% deterministic strategies for Alice and Bob
DA = detPax(nx,na); % size: [na_lambda, na, nx]
DB = detPax(ny,nb); % size: [nb_lambda, nb, ny]

% assemblage reconstruction S{a,b,x,y}
for a = 1:na
  for b = 1:nb
    for x = 1:nx
      for y = 1:ny
        S{a,b,x,y} = 0*sdpvar(d,d);
        for nal = 1:na_lambda
          for nbl = 1:nb_lambda
              i = (nal-1)*nb_lambda + nbl;
              S{a,b,x,y} = S{a,b,x,y} + ...
                   DA(nal,a,x) * DB(nbl,b,y) * s{i};
          end
        end
      end
    end
  end
end

% constraints
sums = 0*sdpvar(d,d);
F = [];

for lambda = 1:nlambda
    F = [F, s{lambda} >= 0];   % σ_λ >= 0
    sums = sums + s{lambda};   % ∑_λ σ_λ
end

for a = 1:na
  for b = 1:nb
    for x = 1:nx
      for y = 1:ny
        F = [F, S{a,b,x,y} - assemblage(:,:,a,b,x,y) >= 0 ];
      end
    end
  end
end

% objective
SR = trace(sums) - 1;

sol = solvesdp(F, SR);

y1 = double(SR);

for lambda = 1:nlambda
    solution_s(:,:,lambda) = double(s{lambda});
end

y2 = solution_s;
disp(sol)

end
