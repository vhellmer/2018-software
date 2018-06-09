% cell_modeling.m
% Rewriting code from Darlington et al., Nat. Comm. 2018
% Simulate metabolic, translational, transcriptional, and orthogonal ribosomal processes 

%% Plot any given variable for each of the four systems
% Key to the variables:
% Access as Y(:, #), where # is the number
% 1 = internalized substrate
% 2 = cell energy
% 3 = cell growth
% 4-7 = mRNAs for transport enzymes, substrate-energy conversion enzymes,  other host proteins, and host ribosomal proteins respectively
% 8 = host rRNA
% 9-12 = loaded ribosomes for transport enzymes, substrate-energy conversion enzymes,  other host proteins, and host ribosomal proteins respectively
% 13-16 = actual proteins for transport enzymes, substrate-energy conversion enzymes,  other host proteins, and host ribosomal proteins respectively
% 17 = host ribosomes
% 18-20 = circuit mRNA, loaded ribosomes, and proteins respectively
% 21 = 16S rRNA (orthogonal rRNA)
% 22 = orthogonal ribosomes


%% Key to the functions:
% initialize("...") sets up all parameters and variables necessary and outputs as [P,Y]. 
% Use "circuit" for including the circuit, "oribo" for including the orthogonal ribosomes and "both" for both
% cell"..."ODE(T,Y,P) returns the derivatives with respect to time for each of the variables as an array dY. 
% Fill in "..." with whatever type of ODE you want to simulate.
%% Quick run how-to:
% To just run the ODE, do the following and replace "..." with the ODE you want to simulate:
% [P,Y] = initialize("...");
% odeoptions = odeset('NonNegative',[1:length(Y)]);
% [T,Y] = ode15s(@(T,Y) cell"..."ODE(T,Y,P), [0 2000], Y, odeoptions);
% Note here that:
% [0 2000] is the time frame you want to simulate over
% T is an array representing time
% Y represents your variables
% P holds your parameters.
% odeoptions ensures that the variables are > 0


clear all; close all; %Clear the screen

%% Shows behavior when new circuit/o-ribosomes introduced into bacteria currently at equilibrium
fplot = figure;
figure(fplot.Number); hold('on');
%set(gca, 'YScale', 'log')
[T,Y] = runAfterSteadyState("circuit"); % run the ODE after reaching steady state with cellplainODE()
for i = 18:19
    plot(T,Y(:,i))
end


function [P,Y] = initialize(options)
    %% initialize(options):
    %   Initializes all the numbers necessary for the model; P contains
    %   all parameters while Y contains all variables
    %   Arguments "options" describes what to add ("circuit" or "ortho" or "both")
    
    %% Constants related to energy
    s_e = 1e4; % external substrate source
    phi_e = 0.5; % efficiency of converting internal substrate into energy
    v_T = 728; % max rate of substrate import
    v_E = 5800; % max rate of substrate-energy conversion
    k_T = 1000; % Michaelis-Menten const. for substrate importer
    k_E = 1000; % Michaelis-Menten const. for substrate-energy conversion

    %% Constants related to transcription
    w_T = 4.14; % max transport enzyme transcription
    w_E = 4.14; % max substrate-energy conversion enzyme transcription
    w_H = 948.93; % max other host protein transcription
    w_R = 930; % max host ribosome transcription
    w_r = 3170; % max host rRNA transcription
    o_T = 4.38; % transcription energy threshold for transporter genes
    o_E = 4.38; % transcription energy threshold for substrate-energy conversion genes
    o_H = 4.38; % transcription energy threshold for other host genes
    o_R = 426.87; % transcription energy threshold for ribosome genes
    o_r = 426.87; % transcription energy threshold for rRNA genes
    k_H = 152219; % other host protein transcription threshold
    h_H = 4; % other host protein transcription Hill constant
    d_m = 0.1; % mRNA degradation rate
    d_r = 0.1; % rRNA degradation rate

    %% Constants related to translation
    b_T = 1; % mRNA-ribosome binding rate
    b_E = 1; % mRNA-ribosome binding rate
    b_H = 1; % mRNA-ribosome binding rate
    b_R = 1; % mRNA-ribosome binding rate
    b_r = 1; % rRNA-ribosome binding rate
    u_T = 1; % mRNA-ribosome unbinding rate
    u_E = 1; % mRNA-ribosome unbinding rate
    u_H = 1; % mRNA-ribosome unbinding rate
    u_R = 1; % mRNA-ribosome unbinding rate
    u_r = 1; % rRNA-ribosome unbinding rate
    d_p = 0; % protein degradation rate
    n_T = 300; % length of transporter enzyme
    n_E = 300; % length of substrate-energy conversion enzyme
    n_H = 300; % length of average other host protein
    n_R = 7459; % length of ribosome
    gamma_max = 1260; % max elongation rate
    k_gamma = 7; % elongation energy threshold
    M_proteome = 1e8; % proteome size
    
    %% Package all parameters together
    P = [s_e;phi_e;v_T;v_E;k_T;k_E;w_T;w_E;w_H;w_R;w_r;o_T;o_E;o_H;o_R;o_r;k_H;h_H;d_m;d_r;b_T;...
        b_E;b_H;b_R;b_r;u_T;u_E;u_H;u_R;u_r;d_p;n_T;n_E;n_H;n_R;gamma_max;k_gamma;M_proteome];
    
    if options == "circuit" % Add circuit constants if necessary
        %% Constants related to circuit genes
        w_Y = 100; % max circuit genes transcription - selected at random/may vary
        o_Y = 4.38; % transcription energy threshold for circuit genes
        b_Y = 1; % mRNA-ribosome binding rate
        u_Y = 1; % mRNA-ribosome binding rate
        n_Y = 300; % length of circuit proteins
        P = [P;w_Y;o_Y;b_Y;u_Y;n_Y];
    end
    
    if options == "oribo" % Add orthogonal ribosome constants if necessary
        %% Constants related to orthogonal ribosome genes
        w_O = 350; % max 16S rRNA genes transcription - selected at random/may vary
        o_O = 4.38; % transcription energy threshold for circuit genes
        b_O = 1; % mRNA-ribosome binding rate
        u_O = 1; % mRNA-ribosome binding rate
        d_O = 0.1; % length of circuit proteins
        P = [P;w_O;o_O;b_O;u_O;d_O];
    end
    
    if options == "both" % Add both
        %% Constants related to circuit genes
        w_Y = 100; % max circuit genes transcription - selected at random/may vary
        o_Y = 4.38; % transcription energy threshold for circuit genes
        b_Y = 1; % mRNA-ribosome binding rate
        u_Y = 1; % mRNA-ribosome binding rate
        n_Y = 300; % length of circuit proteins
        %% Constants related to orthogonal ribosome genes
        w_O = 350; % max 16S rRNA genes transcription - selected at random/may vary
        o_O = 4.38; % transcription energy threshold for circuit genes
        b_O = 1; % mRNA-ribosome binding rate
        u_O = 1; % mRNA-ribosome binding rate
        d_O = 0.1; % length of circuit proteins
        P = [P;w_Y;o_Y;b_Y;u_Y;n_Y;w_O;o_O;b_O;u_O;d_O];
    end
    
    %% Variables related to energy/growth
    s_i = 1000; % internalized substrate
    E_cell = 1000; % cellular energy
    N_cell = 0; % total growth (integral of growth rate)

    %% Variables related to mRNAs
    m_T = 10; % mRNAs for substrate transport enzymes
    m_E = 10; % mRNAs for substrate-energy conversion enzymes
    m_H = 10; % mRNAs for other host proteins
    m_R = 10; % mRNAs for ribosomal proteins
    
    %% Variables related to proteins
    r_R = 10; % rRNAs
    c_T = 10; % translation complexes for substrate transport enzymes
    c_E = 10; % translation complexes for substrate-energy conversion enzymes
    c_H = 10; % translation complexes for other host proteins
    c_R = 10; % translation complexes for ribosomal proteins
    p_T = 10; % substrate transport enzymes
    p_E = 10; % substrate-energy conversion enzymes
    p_H = 10; % other host proteins
    p_R = 10; % ribosomal proteins
    R_h = 10; % functional host ribosomes
   
    %% Package all variables into a single array to use matlab stiff differential equation solver (ode15s)
    Y = [s_i; E_cell; N_cell; m_T; m_E; m_H; m_R; r_R; c_T; c_E; c_H; c_R; p_T; p_E; p_H; p_R; R_h];
    
    if options == "circuit" % Add circuit variables if necessary
        %% Variables for circuit genes
        m_Y = 0; % mRNAs for circuit genes
        c_Y = 0; % translation complexes for circuit genes
        p_Y = 0; % circuit proteins
        Y = [Y;m_Y;c_Y;p_Y];
    end
    if options == "oribo" % Add orthogonal ribosome variables if necessary
        %% Variables for circuit genes
        r_O = 0; % 16S rRNA for orthogonal ribosomes
        R_O = 0; % functional orthogonal ribosomes
        Y = [Y;r_O;R_O];
    end
    if options == "both" % Add circuit variables if necessary
        %% Variables for circuit genes
        m_Y = 0; % mRNAs for circuit genes
        c_Y = 0; % translation complexes for circuit genes
        p_Y = 0; % circuit proteins
        %% Variables for orthogonal ribosome genes
        r_O = 0; % 16S rRNA for orthogonal ribosomes
        R_O = 0; % functional orthogonal ribosomes
        Y = [Y;m_Y;c_Y;p_Y;r_O;R_O];
    end
end


function dY = cellplainODE(T,Y,P)
    %% cellODE(T,Y,P,options):
    %   Calculates derivative with respect to time of all the variable in
    %   the simulation using equations from the Nat. Comm. 2018 paper
    %   No circuit genes or orthogonal ribosomes
    
    %% Unpack parameters from P
    s_e=P(1);phi_e=P(2);v_T=P(3);v_E=P(4);k_T=P(5);k_E=P(6);w_T=P(7);w_E=P(8);
    w_H=P(9);w_R=P(10);w_r=P(11);o_T=P(12);o_E=P(13);o_H=P(14);o_R=P(15);o_r=P(16);
    k_H=P(17);h_H=P(18);d_m=P(19);d_r=P(20);b_T=P(21);b_E=P(22);b_H=P(23);b_R=P(24);
    b_r=P(25);u_T=P(26);u_E=P(27);u_H=P(28);u_R=P(29);u_r=P(30);d_p=P(31);n_T=P(32);n_E=P(33);
    n_H=P(34);n_R=P(35);gamma_max=P(36);k_gamma=P(37);M_proteome=P(38);
    
    %% Unpack current variables from Y
    s_i = Y(1); E_cell = Y(2); N_cell = Y(3); m_T = Y(4); m_E = Y(5); m_H = Y(6);
    m_R = Y(7); r_R = Y(8); c_T = Y(9); c_E = Y(10); c_H = Y(11); c_R = Y(12); 
    p_T = Y(13); p_E = Y(14); p_H = Y(15); p_R = Y(16); R_h = Y(17);

    %% Determine global translation rate (auxiliary variable)
    gamma = gamma_max*E_cell/(k_gamma+E_cell);
    
    %% Determine growth rate (auxiliary variable)
    lambda = gamma*(c_T+c_E+c_H+c_R)/M_proteome; % Supp. eqn 11
    
    %% Determine rate of change for energy/growth variables
    ds_i = v_T*p_T*s_e/(k_T+s_e) - v_E*p_E*s_i/(k_E+s_i) - lambda*s_i; % Supp. eqn 1
    dE_cell = phi_e*v_E*p_E*s_i/(k_E+s_i) - gamma*(c_T+c_E+c_H+c_R) - lambda*E_cell; % Supp. eqn 2, with simplification
    dN_cell = lambda;
    
    %% Determine rate of change for mRNA variables
    dm_T = w_T*E_cell/(E_cell+o_T) + gamma*c_T/n_T - b_T*R_h*m_T + u_T*c_T - (d_m+lambda)*m_T; % Supp. eqn 4, with simplification and sign change on last term
    dm_E = w_E*E_cell/(E_cell+o_E) + gamma*c_E/n_E - b_E*R_h*m_E + u_E*c_E - (d_m+lambda)*m_E; % Same as above
    dm_R = w_R*E_cell/(E_cell+o_R) + gamma*c_R/n_R - b_R*R_h*m_R + u_R*c_R - (d_m+lambda)*m_R; % Same as above
    dm_H = w_H/(1+(p_H/k_H)^h_H)*E_cell/(E_cell+o_H) + gamma*c_H/n_H ...
            - b_H*R_h*m_H + u_H*c_H - (d_m+lambda)*m_H; % Supp. eqn 4, including correct R factor for regulation of host proteins
    
    %% Determine rate of change for proteins
    dc_T = b_T*R_h*m_T - u_T*c_T - gamma*c_T/n_T - (d_p+lambda)*c_T; % Supp. eqn 5, with simplification
    dc_E = b_E*R_h*m_E - u_E*c_E - gamma*c_E/n_E - (d_p+lambda)*c_E; % Same as above
    dc_R = b_R*R_h*m_R - u_R*c_R - gamma*c_R/n_R - (d_p+lambda)*c_R; % Same as above
    dc_H = b_H*R_h*m_H - u_H*c_H - gamma*c_H/n_H - (d_p+lambda)*c_H; % Same as above
    dp_T = gamma*c_T/n_T - (d_p+lambda)*p_T; % Supp. eqn 7, with simplification
    dp_E = gamma*c_E/n_E - (d_p+lambda)*p_E; % Same as above
    dp_H = gamma*c_H/n_H - (d_p+lambda)*p_H; % Same as above
    dp_R = gamma*c_R/n_R - (d_p+lambda)*p_R - b_r*p_R*r_R + u_r*R_h; % Supp. eqn 9
    dr_R = w_r*E_cell/(E_cell+o_r) - b_r*p_R*r_R + u_r*R_h - (d_r+lambda)*r_R; % Supp. eqn 8
    dR_h = b_r*p_R*r_R - u_r*R_h  - (d_p+lambda)*R_h + gamma*c_T/n_T+gamma*c_E/n_E+gamma*c_H/n_H+gamma*c_R/n_R + u_T*c_T+u_E*c_E+u_R*c_R+u_H*c_H ...
            - (b_T*R_h*m_T+b_E*R_h*m_E+b_R*R_h*m_R+b_H*R_h*m_H); % Supp. eqn 10, with sign change on the sum to fix error
    
    dY = [ds_i; dE_cell; dN_cell; dm_T; dm_E; dm_H; dm_R; dr_R; dc_T; dc_E; dc_H; dc_R; dp_T; dp_E; dp_H; dp_R; dR_h];

end


function dY = cellcircuitODE(T,Y,P)
    %% cellcircuitODE(T,Y,P,options):
    %   Calculates derivative with respect to time of all the variable in
    %   the simulation using equations from the Nat. Comm. 2018 paper
    %   No orthogonal ribosomes, but circuit genes included
    
    %% Unpack parameters from P
    s_e=P(1);phi_e=P(2);v_T=P(3);v_E=P(4);k_T=P(5);k_E=P(6);w_T=P(7);w_E=P(8);
    w_H=P(9);w_R=P(10);w_r=P(11);o_T=P(12);o_E=P(13);o_H=P(14);o_R=P(15);o_r=P(16);
    k_H=P(17);h_H=P(18);d_m=P(19);d_r=P(20);b_T=P(21);b_E=P(22);b_H=P(23);b_R=P(24);
    b_r=P(25);u_T=P(26);u_E=P(27);u_H=P(28);u_R=P(29);u_r=P(30);d_p=P(31);n_T=P(32);n_E=P(33);
    n_H=P(34);n_R=P(35);gamma_max=P(36);k_gamma=P(37);M_proteome=P(38);w_Y=P(39);o_Y=P(40);
    b_Y=P(41);u_Y=P(42);n_Y=P(43);
    
    %% Unpack current variables from Y
    s_i = Y(1); E_cell = Y(2); N_cell = Y(3); m_T = Y(4); m_E = Y(5); m_H = Y(6);
    m_R = Y(7); r_R = Y(8); c_T = Y(9); c_E = Y(10); c_H = Y(11); c_R = Y(12); 
    p_T = Y(13); p_E = Y(14); p_H = Y(15); p_R = Y(16); R_h = Y(17);m_Y = Y(18);
    c_Y = Y(19);p_Y = Y(20);

    %% Determine global translation rate (auxiliary variable)
    gamma = gamma_max*E_cell/(k_gamma+E_cell);
    
    %% Determine growth rate (auxiliary variable)
    lambda = gamma*(c_T+c_E+c_H+c_R+c_Y)/M_proteome; % Supp. eqn 11
    
    %% Determine rate of change for energy/growth variables
    ds_i = v_T*p_T*s_e/(k_T+s_e) - v_E*p_E*s_i/(k_E+s_i) - lambda*s_i; % Supp. eqn 1
    dE_cell = phi_e*v_E*p_E*s_i/(k_E+s_i) - gamma*(c_T+c_E+c_H+c_R+c_Y) - lambda*E_cell; % Supp. eqn 12, with simplification
    dN_cell = lambda;
    
    %% Determine rate of change for mRNA variables
    dm_T = w_T*E_cell/(E_cell+o_T) + gamma*c_T/n_T - b_T*R_h*m_T + u_T*c_T - (d_m+lambda)*m_T; % Supp. eqn 4, with simplification and sign change on last term
    dm_E = w_E*E_cell/(E_cell+o_E) + gamma*c_E/n_E - b_E*R_h*m_E + u_E*c_E - (d_m+lambda)*m_E; % Same as above
    dm_R = w_R*E_cell/(E_cell+o_R) + gamma*c_R/n_R - b_R*R_h*m_R + u_R*c_R - (d_m+lambda)*m_R; % Same as above
    dm_Y = w_Y*E_cell/(E_cell+o_Y) + gamma*c_Y/n_Y - b_Y*R_h*m_Y + u_Y*c_Y - (d_m+lambda)*m_Y; % Same as above, but for circuit genes
    dm_H = w_H/(1+(p_H/k_H)^h_H)*E_cell/(E_cell+o_H) + gamma*c_H/n_H - b_H*R_h*m_H + u_H*c_H ...
            - (d_m+lambda)*m_H; % Supp. eqn 4, including correct R factor for regulation of host proteins
    
    %% Determine rate of change for proteins
    dc_T = b_T*R_h*m_T - u_T*c_T - gamma*c_T/n_T - (d_p+lambda)*c_T; % Supp. eqn 5, with simplification
    dc_E = b_E*R_h*m_E - u_E*c_E - gamma*c_E/n_E - (d_p+lambda)*c_E; % Same as above
    dc_R = b_R*R_h*m_R - u_R*c_R - gamma*c_R/n_R - (d_p+lambda)*c_R; % Same as above
    dc_H = b_H*R_h*m_H - u_H*c_H - gamma*c_H/n_H - (d_p+lambda)*c_H; % Same as above
    dc_Y = b_Y*R_h*m_Y - u_Y*c_Y - gamma*c_Y/n_Y - (d_p+lambda)*c_Y; % Same as above, but for circuit genes
    dp_T = gamma*c_T/n_T - (d_p+lambda)*p_T; % Supp. eqn 7, with simplification
    dp_E = gamma*c_E/n_E - (d_p+lambda)*p_E; % Same as above
    dp_H = gamma*c_H/n_H - (d_p+lambda)*p_H; % Same as above
    dp_Y = gamma*c_Y/n_Y - (d_p+lambda)*p_Y; % Same as above, but for circuit genes
    dp_R = gamma*c_R/n_R - (d_p+lambda)*p_R - b_r*p_R*r_R + u_r*R_h; % Supp. eqn 9
    dr_R = w_r*E_cell/(E_cell+o_r) - b_r*p_R*r_R + u_r*R_h - (d_r+lambda)*r_R; % Supp. eqn 8
    dR_h = b_r*p_R*r_R - u_r*R_h  - (d_p+lambda)*R_h + gamma*(c_T/n_T+c_E/n_E+c_H/n_H+c_R/n_R+c_Y/n_Y) + u_T*c_T+u_E*c_E+u_R*c_R+u_H*c_H+u_Y*c_Y ...
            - (b_T*R_h*m_T+b_E*R_h*m_E+b_R*R_h*m_R+b_H*R_h*m_H+b_Y*R_h*m_Y); % Supp. eqn 13, with sign change on the sum to fix error
    
    dY = [ds_i; dE_cell; dN_cell; dm_T; dm_E; dm_H; dm_R; dr_R; dc_T; dc_E; dc_H; dc_R; dp_T; dp_E; dp_H; dp_R; dR_h;dm_Y;dc_Y;dp_Y];

end


function dY = celloriboODE(T,Y,P)
    %% celloriboODE(T,Y,P,options):
    %   Calculates derivative with respect to time of all the variable in
    %   the simulation using equations from the Nat. Comm. 2018 paper
    %   Orthogonal ribosomes included, but these ribosomes are not active
    
    %% Unpack parameters from P
    s_e=P(1);phi_e=P(2);v_T=P(3);v_E=P(4);k_T=P(5);k_E=P(6);w_T=P(7);w_E=P(8);
    w_H=P(9);w_R=P(10);w_r=P(11);o_T=P(12);o_E=P(13);o_H=P(14);o_R=P(15);o_r=P(16);
    k_H=P(17);h_H=P(18);d_m=P(19);d_r=P(20);b_T=P(21);b_E=P(22);b_H=P(23);b_R=P(24);
    b_r=P(25);u_T=P(26);u_E=P(27);u_H=P(28);u_R=P(29);u_r=P(30);d_p=P(31);n_T=P(32);n_E=P(33);
    n_H=P(34);n_R=P(35);gamma_max=P(36);k_gamma=P(37);M_proteome=P(38);w_O=P(39);o_O=P(40);b_O=P(41);u_O=P(42);d_O=P(43);
    
    %% Unpack current variables from Y
    s_i = Y(1); E_cell = Y(2); N_cell = Y(3); m_T = Y(4); m_E = Y(5); m_H = Y(6);
    m_R = Y(7); r_R = Y(8); c_T = Y(9); c_E = Y(10); c_H = Y(11); c_R = Y(12); 
    p_T = Y(13); p_E = Y(14); p_H = Y(15); p_R = Y(16); R_h = Y(17);r_O = Y(18); R_O = Y(19);

    %% Determine global translation rate (auxiliary variable)
    gamma = gamma_max*E_cell/(k_gamma+E_cell);
    
    %% Determine growth rate (auxiliary variable)
    lambda = gamma*(c_T+c_E+c_H+c_R)/M_proteome; % Supp. eqn 11
    
    %% Determine rate of change for energy/growth variables
    ds_i = v_T*p_T*s_e/(k_T+s_e) - v_E*p_E*s_i/(k_E+s_i) - lambda*s_i; % Supp. eqn 1
    dE_cell = phi_e*v_E*p_E*s_i/(k_E+s_i) - gamma*(c_T+c_E+c_H+c_R) - lambda*E_cell; % Supp. eqn 2, with simplification
    dN_cell = lambda;
    
    %% Determine rate of change for mRNA variables
    dm_T = w_T*E_cell/(E_cell+o_T) + gamma*c_T/n_T - b_T*R_h*m_T + u_T*c_T - (d_m+lambda)*m_T; % Supp. eqn 4, with simplification and sign change on last term
    dm_E = w_E*E_cell/(E_cell+o_E) + gamma*c_E/n_E - b_E*R_h*m_E + u_E*c_E - (d_m+lambda)*m_E; % Same as above
    dm_R = w_R*E_cell/(E_cell+o_R) + gamma*c_R/n_R - b_R*R_h*m_R + u_R*c_R - (d_m+lambda)*m_R; % Same as above
    dm_H = w_H/(1+(p_H/k_H)^h_H)*E_cell/(E_cell+o_H) + gamma*c_H/n_H - b_H*R_h*m_H + u_H*c_H ...
            - (d_m+lambda)*m_H; % Supp. eqn 4, including correct R factor for regulation of host proteins
    
    %% Determine rate of change for proteins
    dc_T = b_T*R_h*m_T - u_T*c_T - gamma*c_T/n_T - (d_p+lambda)*c_T; % Supp. eqn 5, with simplification
    dc_E = b_E*R_h*m_E - u_E*c_E - gamma*c_E/n_E - (d_p+lambda)*c_E; % Same as above
    dc_R = b_R*R_h*m_R - u_R*c_R - gamma*c_R/n_R - (d_p+lambda)*c_R; % Same as above
    dc_H = b_H*R_h*m_H - u_H*c_H - gamma*c_H/n_H - (d_p+lambda)*c_H; % Same as above
    dp_T = gamma*c_T/n_T - (d_p+lambda)*p_T; % Supp. eqn 7, with simplification
    dp_E = gamma*c_E/n_E - (d_p+lambda)*p_E; % Same as above
    dp_H = gamma*c_H/n_H - (d_p+lambda)*p_H; % Same as above
    dp_R = gamma*c_R/n_R - (d_p+lambda)*p_R - b_r*p_R*r_R + u_r*R_h - b_O*p_R*r_O + u_O*R_O; % Supp. eqn 16
    dr_R = w_r*E_cell/(E_cell+o_r) - b_r*p_R*r_R + u_r*R_h - (d_r+lambda)*r_R; % Supp. eqn 8
    dr_O = w_O*E_cell/(E_cell+o_O) - b_O*p_R*r_O + u_O*R_O - (d_O+lambda)*r_O; % Supp. eqn 15
    dR_O = b_O*p_R*r_O - u_O*R_O - (d_p+lambda)*R_O; % Supp. eqn 17
    dR_h = b_r*p_R*r_R - u_r*R_h  - (d_p+lambda)*R_h + gamma*c_T/n_T+gamma*c_E/n_E+gamma*c_H/n_H+gamma*c_R/n_R + u_T*c_T+u_E*c_E+u_R*c_R+u_H*c_H ...
            - (b_T*R_h*m_T+b_E*R_h*m_E+b_R*R_h*m_R+b_H*R_h*m_H); % Supp. eqn 10, with sign change on the sum to fix error
    
    dY = [ds_i; dE_cell; dN_cell; dm_T; dm_E; dm_H; dm_R; dr_R; dc_T; dc_E; dc_H; dc_R; dp_T; dp_E; dp_H; dp_R; dR_h; dr_O; dR_O];

end


function dY = cellbothODE(T,Y,P)
    %% cellbothODE(T,Y,P,options):
    %   Calculates derivative with respect to time of all the variable in
    %   the simulation using equations from the Nat. Comm. 2018 paper
    %   Both circuit genes and o-ribosomes included
    %   Circuit genes only translated by o-ribosomes
    
    %% Unpack parameters from P
    s_e=P(1);phi_e=P(2);v_T=P(3);v_E=P(4);k_T=P(5);k_E=P(6);w_T=P(7);w_E=P(8);
    w_H=P(9);w_R=P(10);w_r=P(11);o_T=P(12);o_E=P(13);o_H=P(14);o_R=P(15);o_r=P(16);
    k_H=P(17);h_H=P(18);d_m=P(19);d_r=P(20);b_T=P(21);b_E=P(22);b_H=P(23);b_R=P(24);
    b_r=P(25);u_T=P(26);u_E=P(27);u_H=P(28);u_R=P(29);u_r=P(30);d_p=P(31);n_T=P(32);n_E=P(33);
    n_H=P(34);n_R=P(35);gamma_max=P(36);k_gamma=P(37);M_proteome=P(38);w_Y=P(39);o_Y=P(40);
    b_Y=P(41);u_Y=P(42);n_Y=P(43);w_O=P(44);o_O=P(45);b_O=P(46);u_O=P(47);d_O=P(48);
    
    %% Unpack current variables from Y
    s_i = Y(1); E_cell = Y(2); N_cell = Y(3); m_T = Y(4); m_E = Y(5); m_H = Y(6);
    m_R = Y(7); r_R = Y(8); c_T = Y(9); c_E = Y(10); c_H = Y(11); c_R = Y(12); 
    p_T = Y(13); p_E = Y(14); p_H = Y(15); p_R = Y(16); R_h = Y(17);m_Y = Y(18);
    c_Y = Y(19);p_Y = Y(20);r_O = Y(21); R_O = Y(22);

    %% Determine global translation rate (auxiliary variable)
    gamma = gamma_max*E_cell/(k_gamma+E_cell);
    
    %% Determine growth rate (auxiliary variable)
    lambda = gamma*(c_T+c_E+c_H+c_R+c_Y)/M_proteome; % Supp. eqn 11
    
    %% Determine rate of change for energy/growth variables
    ds_i = v_T*p_T*s_e/(k_T+s_e) - v_E*p_E*s_i/(k_E+s_i) - lambda*s_i; % Supp. eqn 1
    dE_cell = phi_e*v_E*p_E*s_i/(k_E+s_i) - gamma*(c_T+c_E+c_H+c_R+c_Y) - lambda*E_cell; % Supp. eqn 12, with simplification
    dN_cell = lambda;
    
    %% Determine rate of change for mRNA variables
    dm_T = w_T*E_cell/(E_cell+o_T) + gamma*c_T/n_T - b_T*R_h*m_T + u_T*c_T - (d_m+lambda)*m_T; % Supp. eqn 4, with simplification and sign change on last term
    dm_E = w_E*E_cell/(E_cell+o_E) + gamma*c_E/n_E - b_E*R_h*m_E + u_E*c_E - (d_m+lambda)*m_E; % Same as above
    dm_R = w_R*E_cell/(E_cell+o_R) + gamma*c_R/n_R - b_R*R_h*m_R + u_R*c_R - (d_m+lambda)*m_R; % Same as above
    dm_Y = w_Y*E_cell/(E_cell+o_Y) + gamma*c_Y/n_Y - b_Y*R_O*m_Y + u_Y*c_Y - (d_m+lambda)*m_Y; % Same as above, but for circuit genes and using o-ribosomes
    dm_H = w_H/(1+(p_H/k_H)^h_H)*E_cell/(E_cell+o_H) + gamma*c_H/n_H - b_H*R_h*m_H + u_H*c_H ...
            - (d_m+lambda)*m_H; % Supp. eqn 4, including correct R factor for regulation of host proteins
    
    %% Determine rate of change for proteins
    dc_T = b_T*R_h*m_T - u_T*c_T - gamma*c_T/n_T - (d_p+lambda)*c_T; % Supp. eqn 5, with simplification
    dc_E = b_E*R_h*m_E - u_E*c_E - gamma*c_E/n_E - (d_p+lambda)*c_E; % Same as above
    dc_R = b_R*R_h*m_R - u_R*c_R - gamma*c_R/n_R - (d_p+lambda)*c_R; % Same as above
    dc_H = b_H*R_h*m_H - u_H*c_H - gamma*c_H/n_H - (d_p+lambda)*c_H; % Same as above
    dc_Y = b_Y*R_O*m_Y - u_Y*c_Y - gamma*c_Y/n_Y - (d_p+lambda)*c_Y; % Same as above, but for circuit genes and using o-ribosomes
    dp_T = gamma*c_T/n_T - (d_p+lambda)*p_T; % Supp. eqn 7, with simplification
    dp_E = gamma*c_E/n_E - (d_p+lambda)*p_E; % Same as above
    dp_H = gamma*c_H/n_H - (d_p+lambda)*p_H; % Same as above
    dp_Y = gamma*c_Y/n_Y - (d_p+lambda)*p_Y; % Same as above, but for circuit genes
    dp_R = gamma*c_R/n_R - (d_p+lambda)*p_R - b_r*p_R*r_R + u_r*R_h - b_O*p_R*r_O + u_O*R_O; % Supp. eqn 16
    dr_R = w_r*E_cell/(E_cell+o_r) - b_r*p_R*r_R + u_r*R_h - (d_r+lambda)*r_R; % Supp. eqn 8
    dr_O = w_O*E_cell/(E_cell+o_O) - b_O*p_R*r_O + u_O*R_O - (d_O+lambda)*r_O; % Supp. eqn 15
    dR_O = b_O*p_R*r_O - u_O*R_O - (d_p+lambda)*R_O + gamma*c_Y/n_Y ...
            - b_Y*R_O*m_Y + u_Y*c_Y; % Supp. eqn 17, with the (corrected) modifications from supp. eqn 13 to allow translation of circuit genes
    dR_h = b_r*p_R*r_R - u_r*R_h  - (d_p+lambda)*R_h + gamma*(c_T/n_T+c_E/n_E+c_H/n_H+c_R/n_R) + u_T*c_T+u_E*c_E+u_R*c_R+u_H*c_H ...
            - (b_T*R_h*m_T+b_E*R_h*m_E+b_R*R_h*m_R+b_H*R_h*m_H); % Supp. eqn 10, with sign change on the sum to fix error
    
    dY = [ds_i; dE_cell; dN_cell; dm_T; dm_E; dm_H; dm_R; dr_R; dc_T; dc_E; dc_H; dc_R; dp_T; dp_E; dp_H; dp_R; dR_h;dm_Y;dc_Y;dp_Y;dr_O;dR_O];

end


function Yf = runToSteadyState()
    %% runToSteadyState(P,Y)
    % Runs the plain ODE until E_cell stabilizes and output only the final values.
    % Allows us to introduce circuit/o-ribosomes/both after equilibrium has been reached already for the cell.
    [P,Y] = initialize("plain"); % initialize parameters/variables
    threshold = 0.1; % threshold for the max dE_cell/dt allowed
    time = 2000; % move forward by this many minutes each round of the while loop
    odeoptions = odeset('NonNegative',[1:length(Y)]);
    [T,Y] = ode15s(@(T,Y) cellplainODE(T,Y,P), [0 time], Y, odeoptions); % initial run
    test = 1;
    while test > 0
        dY = cellplainODE(T,Y,P);
        if dY(2) < threshold % if dE_cell/dt is smaller than threshold, move forward "time" minutes again and we are done.
            [T,Y] = ode15s(@(T,Y) cellplainODE(T,Y,P), [0 time], Y(length(Y),:), odeoptions);
            test = 0;
        end
        [T,Y] = ode15s(@(T,Y) cellplainODE(T,Y,P), [0 time], Y(length(Y),:), odeoptions); % otherwise, simulate another "time" minutes
    end
    Yf = Y(end,:); % Return final values of variables after steady-state reached
end

function [T,Y] = runAfterSteadyState(options)
    %% Initialize and run cellplainODE() to get steady state in the absence of other genes.
    %% Then, run the ODE described in options ("circuit"/"oribo"/"both") and output the results.
    [P,Y] = initialize(options); % Initialize variables/parameters
    Yf = runToSteadyState(); % Run until steady state and use as initial condition
    switch options
        case "circuit"
            ODE = @cellcircuitODE;
        case "oribo"
            ODE = @celloriboODE;
        case "both"
            ODE = @cellbothODE;
        otherwise
            ODE = @cellplainODE;
    end
    for i = [length(Yf)+1:length(Y)]  % Add any remaining initial variables required for the ODE to Yf
        Yf = [Yf Y(i)];
    end
    odeoptions = odeset('NonNegative',[1:length(Yf)]);
    [T,Y] = ode15s(@(T,Y) ODE(T,Y,P), [0 2000], Yf, odeoptions); % Run the correct ODE
end







