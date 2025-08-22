function [QA_val, Ak_val] = adversary_chk(L, A, w, al, E, delta, zM)
    % INPUT
    % L: in-degree Laplacian matrix
    % A: adjacency matrix
    % w: weighting performance
    % al: number of attack nodes
    % E: maximum attack energy
    % delta: alarm threshold
    % zM: monitor set for CH-k defender
    % OUTPUT
    % QA_val: maximum disruption
    % Ak_val: best attack set
    N = size(L, 1);
    be = length(find(zM == 1)); % number of monitor nodes
    A_adm_set = nchoosek(1:N, al); % all the indices of combination of admissible attack nodes
    nA_adm_set = size(A_adm_set, 1); % number of attack nodes
    Qk = sdpvar(1, 1); % maximum attack impact
    gmA = cell(nA_adm_set, 1);
    psiA = cell(nA_adm_set, 1);
    PA = cell(nA_adm_set, 1);
    eps = cell(nA_adm_set, 1); % epsilon-bound

    % Find the maximum impact and the CH-k policy (at the same time)
    F = [];
    h = 0;
    F = [F, Qk >= 0];
    for nA = 1 : nA_adm_set 
        gmA{nA, 1} = sdpvar(N, 1);
        psiA{nA, 1} = sdpvar(al, 1);
        eps{nA, 1} = sdpvar(1, 1);
        PA{nA, 1} = sdpvar(N, N, 'symmetric');       
        BA = zeros(N, al);
        for i = 1 : al
            BA(:, i) = A(:, A_adm_set(nA, i));
        end
        lmi = [-L'*PA{nA, 1}-PA{nA, 1}*L+diag(w) PA{nA, 1}*BA;
               BA'*PA{nA, 1} -diag(psiA{nA, 1})] - diag([gmA{nA, 1}.*zM; zeros(al, 1)]);
        F = [F, lmi <= 0];
        F = [F, gmA{nA, 1} >= 0];
        F = [F, psiA{nA, 1} >= 0];
        F = [F, eps{nA, 1} >= 0];
        F = [F, eps{nA, 1} <= 1e-1];
        F = [F, delta'*gmA{nA, 1} + E * ones(1, al) * psiA{nA, 1} <= Qk - eps{nA, 1}];
        h = h - eps{nA, 1};
    end
    h = h + 1e3*Qk;
    options = sdpsettings('verbose',0,'debug',0,'solver','mosek');
    sol = optimize(F, h, options);
    % disp(sol.solvertime);
    QA_val = value(Qk);

    % Find the CH-k adversary policy Ak
    % eps = cell(nA_adm_set, 1); % epsilon-bound
    % F = [];
    % h = 0;
    % for nA = 1 : nA_adm_set
    %     eps{nA, 1} = sdpvar(1, 1);
    %     BA = zeros(N, al);
    %     for i = 1 : al
    %         BA(:, i) = A(:, A_adm_set(nA, i));
    %     end
    %     lmi = [-L'*PA{nA, 1}-PA{nA, 1}*L+diag(w) PA{nA, 1}*BA;
    %            BA'*PA{nA, 1} -diag(psiA{nA, 1})] - diag([gmA{nA, 1}.*zM; zeros(al, 1)]);
    %     F = [F, lmi <= 0];
    %     F = [F, gmA{nA, 1} >= 0];
    %     F = [F, psiA{nA, 1} >= 0];
    %     F = [F, eps{nA, 1} >= 0];
    %     F = [F, delta'*gmA{nA, 1} + E * ones(1, al) * psiA{nA, 1} <= Q_val - eps{nA, 1}];
    %     h = h - eps{nA, 1};
    % end
    % options = sdpsettings('verbose',2,'debug',0,'solver','mosek');
    % optimize(F, h, options);
    for nA = 1 : nA_adm_set  
         % disp(value(eps{nA, 1}));
        if value(eps{nA, 1}) <= 1e-3
            Ak_val = A_adm_set(nA, :);
            %break;
        end
    end
end