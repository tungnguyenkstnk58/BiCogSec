function [QM_val, zM_val] = defender_chk(L, A, w, al, E, delta, Ak, pk, cs, be, chk)
    % INPUT
    % L: in-degree Laplacian matrix
    % A: adjacency matrix
    % w: weighting performance
    % al: number of attack nodes
    % E: maximum attack energy
    % delta: alarm threshold
    % Ak: CH-k adversary policies
    % pk: belief in CH-k
    % cs: sensor cost
    % be: maxmimum monitor nodes
    % chk: CH-k defender policy
    % OUTPUT
    % QM_val: minimum disruption
    % zM_val: best monitor set
    N = size(L, 1);
    A_adm_set = nchoosek(1:N, al); % all the indices of combination of admissible attack nodes
    nA_adm_set = size(A_adm_set, 1); % number of attack nodes
    zM = binvar(N, 1); % monitor set
    zM_val = zeros(N, 1);
    
    % Find the CH-1 defender policy
    if chk == 1
        Qk = sdpvar(1, 1); % maximum attack impact
        gmA = cell(nA_adm_set, 1);
        psiA = cell(nA_adm_set, 1);
        PA = cell(nA_adm_set, 1);
        F = [];
        F = [F, Qk >= 0];
        F = [F, ones(1, N) * zM <= be];
        for nA = 1 : nA_adm_set 
            gmA{nA, 1} = sdpvar(N, 1);
            psiA{nA, 1} = sdpvar(al, 1);
            PA{nA, 1} = sdpvar(N, N, 'symmetric');       
            BA = zeros(N, al);
            for i = 1 : al
                BA(:, i) = A(:, A_adm_set(nA, i));
            end
            lmi = [-L'*PA{nA, 1}-PA{nA, 1}*L+diag(w) PA{nA, 1}*BA;
                    BA'*PA{nA, 1} -diag(psiA{nA, 1})] - diag([gmA{nA, 1}; zeros(al, 1)]);
            F = [F, lmi <= 0];
            F = [F, gmA{nA, 1} >= 0];
            F = [F, gmA{nA, 1} <= 1e5 * zM];
            F = [F, psiA{nA, 1} >= 0];
            F = [F, delta'*gmA{nA, 1} + E * ones(1, al) * psiA{nA, 1} <= Qk];
        end
        h = cs' * zM + Qk;
        options = sdpsettings('verbose',0,'debug',0,'solver','bnb','bnb.solver','mosek');
        optimize(F, h, options);
        zM_val = value(zM);
        QM_val = value(h);
    end
    % Find the CH-k (k > 1) defender policy
    if chk > 1
        gmA = cell(chk-1, 1);
        psiA = cell(chk-1, 1);
        PA = cell(chk-1, 1);
        F = [];
        % h = 0;
        h = sdpvar(1, 1);
        F = [F, ones(1, N) * zM <= be];    
        for chi = 1 : chk-1
            gmA{chi, 1} = sdpvar(N, 1);
            psiA{chi, 1} = sdpvar(al, 1);
            PA{chi, 1} = sdpvar(N, N, 'symmetric');
            BA = zeros(N, al);
            for i = 1 : al
                BA(:, i) = A(:, Ak{chi, 1}(i));
            end
            lmi = [-L'*PA{chi, 1}-PA{chi, 1}*L+diag(w) PA{chi, 1}*BA;
                    BA'*PA{chi, 1} -diag(psiA{chi, 1})] - diag([gmA{chi, 1}; zeros(al, 1)]);
            F = [F, lmi <= 0];
            F = [F, gmA{chi, 1} >= 0];
            F = [F, gmA{chi, 1} <= 1e5 * zM];
            F = [F, psiA{chi, 1} >= 0];
            F = [F, delta'*gmA{chi, 1} + E * ones(1, al) * psiA{chi, 1} <= h];
            % h = h + pk(chi) * (delta'*gmA{chi, 1} + E * ones(1, al) * psiA{chi, 1});
        end
        obj = cs' * zM + h;
        options = sdpsettings('verbose',0,'debug',0,'solver','bnb','bnb.solver','mosek');
        optimize(F, obj, options);
        zM_val = value(zM);
        QM_val = value(obj);
    end
end