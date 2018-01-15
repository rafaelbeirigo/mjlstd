initial_mode = 1;

x0 = [1; 1];

epsilon_F = 1e-9;
epsilon_Y = epsilon_F;
epsilon_sum = epsilon_F;

%% TODO: inicializar adequadamente
gamma_var = 0.1;
lambda_var = 0.1;

%% Number of repetitions of the experiment
num_repetitions = 1e3;

%% Maximum number of steps to simulate
T = 1e3;
