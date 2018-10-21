classdef PhysicsConstants
    %PhysicsConstants.m provides commonly-used physics constants
    % Constants include:
    %
    %
    %  hbar_J     [J*s]            reduced Planck constant
    %
    %  hbar_ev    [eV*s]           reduced Planck constant
    %
    %  c          [m/s]            speed of light
    %
    %  c_nm       [nm/s]           speed of light
    %
    %  q          [C]              elementary charge (this is positive!)
    %
    %  kBeV       [eV/K]           Boltzmann constant in eV
    %
    %  kBJ        [J/K]            Boltzmann constant in eV
    %
    %  eps0       [F/m = C/(V*m)]  Permittivity of free space
    %
    %  r_Bohr_nm  [nm]             Bohr radius
    %
    
    properties (Constant)
        hbar_J = 1.054571726E-34; % [J*s]  reduced Planck constant
        hbar_ev = 6.58211928E-16; % [eV*s] reduced Planck constant
        
        c = 2.99792E8; % [m/s] speed of light
        c_nm = 2.99792E17; % [nm/s] speed of light
        
        q = 1.60217657E-19; % [C] elementary charge
        
        kBeV = 8.6173324E-5; % [eV/K] Boltzmann constant in eV
        kBJ = 1.3806488E-23; % [J/K] Boltzmann constant in eV
        
        eps0 = 8.854E-12; % [F/m = C/(V*m)] Permittivity of free space
        
        r_Bohr_nm = 0.0529177; % [nm] Bohr radius
        
    end
    
end

