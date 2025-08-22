clear
close all
clc
%% System parameters
N = 10;
[L, A] = ConnectedGraph(N, 0); % Laplacian/Adjacency matrix
Theta = 0.5 * eye(N); % self-loop gains
L = L + Theta; % weighted Laplacian with self-loops
delta = 1*ones(N, 1) + 0.1 * rand(N, 1); % alarm threshold
E = 1e3; % Attack maximum energy
al = 2; % attack nodes
be = 2; % monitor nodes
chq = 10; % maximum CH
lamP = 6; % Poisson distribution
ch1_q = 1 : chq; % CH levels
pkr = lamP.^ch1_q * exp(-lamP) ./ factorial(ch1_q); % Poisson distribution for the reference belief
w = ones(N,1) + 0.1 * randn(N, 1); % performance weighting factor
cs = 1*ones(N,1) + 0.1 * randn(N, 1); % cost of sensors

%% Players action
Ak = cell(chq, 1); % Adversary policy CH-1 -> Ch-q
QAk = zeros(chq, 1); % Adversary payoff
Ak_time = zeros(chq, 1); % computation time for Ak
Mk = cell(chq, 1); % Defender policy CH-1 -> CH-q
zMk = cell(chq, 1); % Defender policy CH-1 -> CH-q
QMk = zeros(chq, 1); % Defende payoff
Mk_time = zeros(chq, 1); % computation time for Mk

tic;
[QM_val_ch1, zM_val_ch1] = defender_chk(L, A, w, al, E, delta, [], [], cs, be, 1);
QMk(1, 1) = QM_val_ch1;
Mk_time(1) = toc;
zMk{1, 1} = zM_val_ch1;
Mk_opt = find(zM_val_ch1 == 1)';
zMk{1, 1} = zeros(N, 1); zMk{1, 1}(1:be) = 1;
Mk{1, 1} = find(zMk{1, 1} == 1)';


for chk = 1 : chq
    tic;
    [Q_val, Ak_val] = adversary_chk(L, A, w, al, E, delta, zMk{chk, 1});
    Ak_time(chk) = toc;
    Ak{chk, 1} = Ak_val;
    QAk(chk, 1) = Q_val;
    if chk < chq
        pk = pkr(1 : chk);
        pk = pk ./ sum(pk);
        tic;
        [QM_val_chk, zM_val_chk] = defender_chk(L, A, w, al, E, delta, Ak, pk, cs, be, chk + 1);
        Mk_time(chk + 1) = toc;
        zMk{chk + 1, 1} = zM_val_chk;
        Mk{chk + 1, 1} = find(zMk{chk + 1, 1} == 1)';
        QMk(chk + 1, 1) = QM_val_chk;
    end
    % if Mk{chk, 1} == Mk{chk + 1, 1}
    %     chq = chk;
    %     break;
    % end
end


%% Result display
result = [];
for chk = 1 : chq
    result = [result; [chk, Ak{chk, 1}, QAk(chk, 1), Mk{chk, 1}, QMk(chk, 1)]];
end
disp(round(result));
figure(1)
plot(1 : chq, QAk, 'LineWidth', 2);
hold on
plot(1 : chq, QMk, 'LineWidth', 2);
grid;
legend('Adversary payoff','Defender payoff');
xlabel('CH-k')
print(gcf,'figure/10node_2v2','-djpeg','-r300');


% figure(2)
% stem(1 : chq, pkr, 'LineWidth', 2);
% grid;


