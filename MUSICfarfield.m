close all
clear
clc
colors = {'b', 'r', [0.3, 0.3, 0.3]}; % azul, roxo, preto
%=========================================================================
antennas = [8,32,64]; % número de antenas no array
snapshots = 200; % número de amostras
lambda = 0.02; % comprimento de onda (3.10^8*freq)
d_antenna = lambda/2; % espaçamento entre antenas em comprimentos de onda.
AoA = [-53, -12, 48, 55]; % ângulos de chegada dos alvos (em graus)
users = length(AoA); % número de sinais
d_user = [50, 24, 130, 192]; %distância dos usuários
%=========================================================================

figure; % criar nova figura para cada SNR
hold on; % manter as curvas no mesmo gráfico
for ii = 1:length(antennas)
    M = antennas(ii); % número de antenas

    % Geração de sinais recebidos com SNR específico
    [Y, X, A] = received_signals(M, snapshots, d_antenna, lambda, AoA, users, d_user); % recebidos

    % Cálculo da matriz de covariância média
    Matriz_cov = (Y * Y') / snapshots;

    % Decomposição em autovalores e autovetores
    [eigenvectors, eigenvalues] = eig(Matriz_cov);

    % Subespaço de ruído    
    Un = n_space(eigenvectors, eigenvalues, users); % noise subspace

    % Geração do Pseudo Espectro MUSIC
    theta = -90:1:90; % eixo x
    Pmusic = zeros(size(theta)); % eixo y
    for i = 1:length(theta)
        a = stevec(M, d_antenna, lambda, theta(i));
        Pmusic(i) = 1 / (a' * (Un * Un') * a);
    end

    % Normalização do pseudoespectro
    Pmusic = abs(Pmusic) / max(abs(Pmusic));

    % Plot do Pseudo Espectro MUSIC para cada valor de M
    plot(theta, 10 * log10(Pmusic), 'Color', colors{ii}, ...
        'DisplayName', ['M = ' num2str(M)]);
end

% Adicionar uma entrada adicional para AoA e distância na legenda, sem traço
legend_entry = sprintf('AoA: %s \nDistance: %s (m)', strjoin(string(AoA), ', '), strjoin(string(d_user), ', '));
plot(nan, nan, 'Color', 'none', 'DisplayName', legend_entry); % usa uma linha invisível sem traço


% Configurações do gráfico
xlabel('Angle (degrees)');
ylabel('Pseudo Spectrum (dB)');
title('MUSIC Spectrum (dB)');
grid on;
legend show; % exibe todas as legendas
hold off; % finaliza o gráfico atual


% Funções locais
function steering_vector = stevec(M,d_antenna,lambda,theta)
    gamma = (2*pi * d_antenna)/lambda;
    steering_vector = exp(-1i*gamma*(0:M-1)'*sind(theta));
end

function noise_subspace = n_space(eigenvectors, eigenvalues, users)
    [~, i] = sort(diag(eigenvalues), 'descend');
    eigenvectors = eigenvectors(:, i);
    noise_subspace = eigenvectors(:, users+1:end);
end

function [Y,X,A] = received_signals(M, snapshots, d_antenna, lambda, AoA, users,d_user)

% Definir parâmetros de simulação
P_tx = 0.1; % Potência transmitida em watts
P_noise_dBm = -90; % Potência de ruído em dBm
P_noise = 10^((P_noise_dBm - 30) / 10); % Converte para watts

% Inicializar matriz de canal e pathloss
H = zeros(M, users); % M: número de antenas, users: número de usuários
pathloss = zeros(users, 1);

for s = 1:users
    % Cálculo do pathloss baseado na distância da fonte s ao arranjo
    pathloss(s) = 1 / (d_user(s) ^ 3);

    % Vetor de steering para a fonte s
    A = stevec(M, d_antenna, lambda, AoA(s)); % resposta do array ao AoA

    % Coeficiente de canal para a fonte s
    H(:, s) = pathloss(s) * A;
end

% Gerar sinais transmitidos
X = sqrt(P_tx) * randn(users, snapshots)/ sqrt(2); % Sinais transmitidos
Y = H * X; % Sinais recebidos pelo array

% Gerar ruído branco com a potência ajustada para -110 dBm
Z = sqrt(P_noise) * (randn(M, snapshots) + 1j * randn(M, snapshots)) / sqrt(2); 

% Adicionar ruído ao sinal recebido
Y = Y + Z;

% Cálculo da SNR para verificação
SNR = 10 * log10(sum(abs(Y).^2, 'all') / sum(abs(Z).^2, 'all'));
disp(['SNR = ', num2str(SNR), ' dB']);


end