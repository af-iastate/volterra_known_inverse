clc; close all; clearvars;

LpSp_fileName = getMostRecentSimulation('lpspn5');
BpSp_fileName = getMostRecentSimulation('bpspn5');
LpDp_fileName = getMostRecentSimulation('lpdpn5');


%% Load Data
result = load(LpSp_Filename, 'h_true', ...
    'h_pseudo_g', 'h_rls_offline_g', 'h_rls_online_g',  ...
    'h_pseudo_mln', 'h_rls_offline_mln', 'h_rls_online_mln', ...
    'h_pseudo_two_tone', 'h_rls_offline_two_tone', 'h_rls_online_two_tone' ...
);

LpSp_h_true                     = result.h_true;

LpSp_h_pseudo_g                 = result.h_pseudo_g{6};
LpSp_h_rls_offline_g            = result.h_rls_offline_g{6};
LpSp_h_rls_online_g             = result.h_rls_online_g{6};

LpSp_h_pseudo_mln               = result.h_pseudo_mln{6};
LpSp_h_rls_offline_mln          = result.h_rls_offline_mln{6};
LpSp_h_rls_online_mln           = result.h_rls_online_mln{6};

LpSp_h_pseudo_two_tone          = result.h_pseudo_two_tone{6};
LpSp_h_rls_offline_two_tone     = result.h_rls_offline_two_tone{6};
LpSp_h_rls_online_two_tone      = result.h_rls_online_two_tone{6};

%% --------------------------------------------------------------------
result = load(BpSp_Filename, 'h_true', ...
    'h_pseudo_g', 'h_rls_offline_g', 'h_rls_online_g',  ...
    'h_pseudo_mln', 'h_rls_offline_mln', 'h_rls_online_mln', ...
    'h_pseudo_two_tone', 'h_rls_offline_two_tone', 'h_rls_online_two_tone' ...
);

BpSp_h_true                     = result.h_true;

BpSp_h_pseudo_g                 = result.h_pseudo_g{6};
BpSp_h_rls_offline_g            = result.h_rls_offline_g{6};
BpSp_h_rls_online_g             = result.h_rls_online_g{6};

BpSp_h_pseudo_mln               = result.h_pseudo_mln{6};
BpSp_h_rls_offline_mln          = result.h_rls_offline_mln{6};
BpSp_h_rls_online_mln           = result.h_rls_online_mln{6};

BpSp_h_pseudo_two_tone          = result.h_pseudo_two_tone{6};
BpSp_h_rls_offline_two_tone     = result.h_rls_offline_two_tone{6};
BpSp_h_rls_online_two_tone      = result.h_rls_online_two_tone{6};

% --------------------------------------------------------------------
result = load(LpDp_Filename, 'h_true', ...
    'h_pseudo_g', 'h_rls_offline_g', 'h_rls_online_g',  ...
    'h_pseudo_mln', 'h_rls_offline_mln', 'h_rls_online_mln', ...
    'h_pseudo_two_tone', 'h_rls_offline_two_tone', 'h_rls_online_two_tone' ...
);

LpDp_h_true                     = result.h_true;

LpDp_h_pseudo_g                 = result.h_pseudo_g{6};
LpDp_h_rls_offline_g            = result.h_rls_offline_g{6};
LpDp_h_rls_online_g             = result.h_rls_online_g{6};

LpDp_h_pseudo_mln               = result.h_pseudo_mln{6};
LpDp_h_rls_offline_mln          = result.h_rls_offline_mln{6};
LpDp_h_rls_online_mln           = result.h_rls_online_mln{6};

LpDp_h_pseudo_two_tone          = result.h_pseudo_two_tone{6};
LpDp_h_rls_offline_two_tone     = result.h_rls_offline_two_tone{6};
LpDp_h_rls_online_two_tone      = result.h_rls_online_two_tone{6};


%% Analyze
LpSp_err_pseudo_g               = norm(LpSp_h_pseudo_g - LpSp_h_true);
LpSp_err_rls_offline_g          = norm(LpSp_h_rls_offline_g - LpSp_h_true);
LpSp_err_rls_online_g           = norm(LpSp_h_rls_online_g - LpSp_h_true);

LpSp_err_pseudo_mln             = norm(LpSp_h_pseudo_mln - LpSp_h_true);
LpSp_err_rls_offline_mln        = norm(LpSp_h_rls_offline_mln - LpSp_h_true);
LpSp_err_rls_online_mln         = norm(LpSp_h_rls_online_mln - LpSp_h_true);

LpSp_err_pseudo_two_tone        = norm(LpSp_h_pseudo_two_tone - LpSp_h_true);
LpSp_err_rls_offline_two_tone   = norm(LpSp_h_rls_offline_two_tone - LpSp_h_true);
LpSp_err_rls_online_two_tone    = norm(LpSp_h_rls_online_two_tone - LpSp_h_true);


BpSp_err_pseudo_g               = norm(BpSp_h_pseudo_g - BpSp_h_true);
BpSp_err_rls_offline_g          = norm(BpSp_h_rls_offline_g - BpSp_h_true);
BpSp_err_rls_online_g           = norm(BpSp_h_rls_online_g - BpSp_h_true);

BpSp_err_pseudo_mln             = norm(BpSp_h_pseudo_mln - BpSp_h_true);
BpSp_err_rls_offline_mln        = norm(BpSp_h_rls_offline_mln - BpSp_h_true);
BpSp_err_rls_online_mln         = norm(BpSp_h_rls_online_mln - BpSp_h_true);

BpSp_err_pseudo_two_tone        = norm(BpSp_h_pseudo_two_tone - BpSp_h_true);
BpSp_err_rls_offline_two_tone   = norm(BpSp_h_rls_offline_two_tone - BpSp_h_true);
BpSp_err_rls_online_two_tone    = norm(BpSp_h_rls_online_two_tone - BpSp_h_true);


LpDp_err_pseudo_g               = norm(LpDp_h_pseudo_g - LpDp_h_true);
LpDp_err_rls_offline_g          = norm(LpDp_h_rls_offline_g - LpDp_h_true);
LpDp_err_rls_online_g           = norm(LpDp_h_rls_online_g - LpDp_h_true);

LpDp_err_pseudo_mln             = norm(LpDp_h_pseudo_mln - LpDp_h_true);
LpDp_err_rls_offline_mln        = norm(LpDp_h_rls_offline_mln - LpDp_h_true);
LpDp_err_rls_online_mln         = norm(LpDp_h_rls_online_mln - LpDp_h_true);

LpDp_err_pseudo_two_tone        = norm(LpDp_h_pseudo_two_tone - LpDp_h_true);
LpDp_err_rls_offline_two_tone   = norm(LpDp_h_rls_offline_two_tone - LpDp_h_true);
LpDp_err_rls_online_two_tone    = norm(LpDp_h_rls_online_two_tone - LpDp_h_true);

% -------------------------------------------------------------------------------

db = @(x) 20*log10(x);

tblMatrix = db([
    LpSp_err_pseudo_g, LpSp_err_rls_offline_g, LpSp_err_rls_online_g, ...
    BpSp_err_pseudo_g, BpSp_err_rls_offline_g, BpSp_err_rls_online_g, ...
    LpDp_err_pseudo_g, LpDp_err_rls_offline_g, LpDp_err_rls_online_g;
    ...
    LpSp_err_pseudo_mln, LpSp_err_rls_offline_mln, LpSp_err_rls_online_mln, ...
    0, 0, 0, ...
    LpDp_err_pseudo_mln, LpDp_err_rls_offline_mln, LpDp_err_rls_online_mln;
    ...
    LpSp_err_pseudo_two_tone, LpSp_err_rls_offline_two_tone, LpSp_err_rls_online_two_tone, ...
    BpSp_err_pseudo_two_tone, BpSp_err_rls_offline_two_tone, BpSp_err_rls_online_two_tone, ...
    LpDp_err_pseudo_two_tone, LpDp_err_rls_offline_two_tone, LpDp_err_rls_online_two_tone
]);

tblCells = arrayfun(@(x) num2str(x, '%.2f'), tblMatrix, 'UniformOutput', 0);
tblCells(2,4:6) = {'--', '--', '--'};
[M,N] = size(tblCells);
for ii=1:M
    row = tblCells(ii,:);
    result = join(row, ' & ');
    result = [result{:} ' \\'];
    disp(result)
end

