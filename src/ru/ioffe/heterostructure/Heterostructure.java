package ru.ioffe.heterostructure;

import ru.ioffe.semiconductor.SemiconductorW;

/**
 * @author Evgeny Shevchenko
 */
public class Heterostructure {

    /**
     * @param args the command line arguments
     */
    @SuppressWarnings("unused")
    public static void main(String[] args) {

        // [1] Jmerik et al. Semiconductors Vol. 42 No. 12 (2008)
        // [2] Caro, Schulz, Healy, O'Reilly, J. Appl. Phys. 109, 084110 (2011)
        // [3] Figge, Kroncke, Hommel, Epelbaum, Appl. Phys. Lett. 94, 101915 (2009)
        // [4] Vurgaftman, Meyer, Ram-Mohan, J. Appl. Phys. 89, 11, 5815 (2001)
        // [5] http://www.ioffe.ru/SVA/NSM/Semicond/
        // GaN = SemiconductorW(
        // nrg_gap = 3.42, # Energy gap (eV) [1]
        // alpha = 0.94e-3, # Coefficient for Varshni equation [ ]
        // beta = 791, # Coefficient for Varshni equation [ ]
        // nrg_gap_0 = 3.4975, # Energy gap at T = 0 tailored from [1]
        // diel_st = 10.4, # Electric constant [ ]
        // a = 3.189, # Lattice parameters (angstrom) [2]
        // c = 5.185, # [2]
        // a0 = 3.1872, # Lattice parameters at T = 0 tailored from [2]
        // c0 = 5.1825, # tailored from [2]
        // a_Tch = 868, # Characteristic temperature, K [3]
        // c_Tch = 898, # [3]
        // a_TEC = 6.24e-6, # Thermal expansion coeffcient, K**-1 [3]
        // c_TEC = 5.73e-6, # [3]
        // mass_elx = 0.2, # Efficient electron mass in x direction [5]
        // mass_elz = 0.2, # Efficient electron mass in z direction [5]
        // mass_hhx = 1.6, # Efficient heavy hole mass in x direction [5]
        // mass_hhz = 1.75, # Efficient heavy hole mass in z direction [5]
        // mass_lhx = 1.1, # Efficient light hole mass in x direction [5]
        // mass_lhz = 0.15, # Efficient light hole mass in z direction [5]
        // c11 = 390., # Elastic constant c11 (GPa) [4]
        // c12 = 145., # Elastic constant c12 (GPa) [4]
        // c13 = 106., # Elastic constant c13 (GPa) [4]
        // c33 = 398., # Elastic constant c33 (GPa) [4]
        // c44 = 105., # Elastic constant c44 (GPa) [4]
        // pz13 = -0.45, # Piezoelectric constant pz13 (C/m**2) [2]
        // pz33 = 0.83, # Piezoelectric constant pz33 (C/m**2) [2]
        // pz15 = -0.38, # Piezoelectric constant pz15 (C/m**2) [2]
        // polar_sp = -0.029, # Spontaneous polarization (C/m**2) [4]

        SemiconductorW GaN = new SemiconductorW(3.189, 5.185, 3.1872, 5.1825, // Lattice parameters a, c, a_0, c_0
                868, 898, 6.24e-6, 5.73e-6, // Lattice coefficietnts for a_Tch, c_Tch, a_TEC, c_TEC
                3.42, 3.4975, // energy gap, energy gap at 0 K
                0.94e-3, 791, // alpha, beta for Varshni equation
                0.2, 0.2, // Electron masses x, z
                1.6, 1.75,// Heavy holes masses x, z
                1.1, 0.15,// Light holes masses x, z
                390, 145, 106, 398, 105, // Elastic constants c11, c12, c13, c33, c44
                -0.45, 0.83, -0.38, // Piezoelectric tensor components pz13, pz33, pz15
                -0.029, 10.4, // Spontaneous polarization, dielectric constant
                300 /* Temperature */);

        // AlN = SemiconductorW(
        // nrg_gap = 6.08, # Energy gap (eV) [Jmerik]
        // alpha = 2.63e-3, # coefficient for Varshni equation
        // beta = 2082, # coefficient for Varshni equation
        // nrg_gap_0 = 6.1794, # Energy gap at T = 0
        // diel_st = 8.5, # Electric constant [1]
        // a = 3.11, # Lattice par a at T = 300 K (Angstrom) [1]
        // c = 4.89, # Lattice par c at T = 300 K (Angstrom) [1]
        // a0 = 3.1112, # Lattice par a at T = 0
        // c0 = 4.9807, # Lattice par a at T = 0
        // a_Tch = 1455, # Characteristic temperature, K
        // c_Tch = 1317,
        // a_TEC = 7.1e-6, # Infinity coeffcient, K**-1
        // c_TEC = 5.8e-6,
        // mass_elx = 0.3, # Efficient electron mass in x direction [1]
        // mass_elz = 0.3, # Efficient electron mass in z direction [1]
        // mass_hhx = 10.42, # Efficient heavy hole mass in x direction
        // mass_hhz = 3.53, # Efficient heavy hole mass in z direction [1]
        // mass_lhx = 0.24, # Efficient light hole mass in x direction
        // mass_lhz = 3.53, # Efficient light hole mass in z direction [1]
        // c11 = 410., # Elastic constant c11, GPa
        // c12 = 149., # Elastic constant c12, GPa
        // c13 = 99., # Elastic constant c13, GPa
        // c33 = 389., # Elastic constant c33, GPa
        // c44 = 125., # Elastic constant c44, GPa
        // pz13 = -0.48, # Piezoelectric constant pz13 (C / m**2) [1]
        // pz33 = 1.55, # Piezoelectric constant pz33 (C / m**2) [1]
        // pz15 = 1, # FIXME Piezoelectric constant pz15 (C / m**2) [1]
        // polar_sp = -0.081, # Spontaneous polarization (C / m**2) [1]
        // name = 'AlN')
        SemiconductorW AlN = new SemiconductorW(3.11, 4.89, 3.1112, 4.9807, // Lattice parameters a, c, a_0, c_0
                1455, 1317, 7.1e-6, 5.8e-6, // Lattice coefficietnts for a_Tch, c_Tch, a_TEC, c_TEC
                6.08, 6.1794, // energy gap, energy gap at 0 K
                2.63e-3, 2082, // alpha, beta for Varshni equation
                0.3, 0.3, // Electron masses x, z
                10.42, 3.53,// Heavy holes masses x, z
                0.24, 3.53,// Light holes masses x, z
                410, 149, 99, 398, 125, // Elastic constants c11, c12, c13, c33, c44
                -0.48, 1.55, -1, // Piezoelectric tensor components pz13, pz33, pz15
                -0.081, 8.5, // Spontaneous polarization, dielectric constant
                300 /* Temperature */);

        // CompoundW AlGaN = new CompoundW(0.3, 300, AlN, GaN, 1.1);
    }
}
