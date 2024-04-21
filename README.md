# trilayer_stacked : Tight Binding Calculation of DOS and MHC-DOS kinetic rates in Stacked Trilayer Graphene 

Calculation density of states for ABA/ABC stacked trilayer graphene using a low-energy momentum space model (or Fourier transformed tight-binding model). The DOS files are used as input to the MHC-DOS model for charge transfer rates with redox couple (Ruthenium Hexamine). For rate model see the Julia-based [ElectrochemicalKinetics.jl](https://github.com/BattModels/ElectrochemicalKinetics.jl) package. 

## Contact 

Dr. Stephen Carr : stephen_carr1@brown.edu

Mohammad Babar : mdbabar@umich.edu 

## Code descriptions

Main file descriptions are as follows, other files are either output or supporting function files. 

1. `trilayer_stacking_band_calc.m` generates DOS files for ABA or ABC stacking using tight binding. named `ABC_dos.mat` or `ABA_dos.mat`. Description of input arguments can be found at the beginning of the file.

2. `script.jl` main Julia calculation script that outputs `.mat` file with oxidation (k_ox) and reduction rates (k_red) for specified parameters $A$ , $\lambda$ and $\eta$.

where $\lambda$ = reorganization energy (eV), $\eta$ = applied overpotential (V) and $A$ = proportionality constant for MHC-DOS theory. Other input parameters:

i. `C_dl` : EDL capacitance (F)

ii. `V_dl` : EDL voltage (V)

iii. `C_q` : Quantum capacitance (F)

iv. `V_q` : Quantum capacitance voltage (V)

v. `Vq_min / Vq_max` : Min/Max range of Quantum capacitance voltage for interpolation (Eq. 3 in paper)

vi. `kT` : Thermal energy to temperature setting (0.26 eV at 300 K)

vii. `ef` : Fermi energy of the electrode

Run Julia scripts using:

```
> julia script.jl
```

The output prints the rates, e.g.

```
ABA k_ox 2.6323040539732973e-5 k_red 2.6323040539732787e-5
ABC k_ox 2.77951068679322e-5 k_red 2.7795106867932337e-5
```

Here LHS and RHS rates are equal as `η = 0.0` (no driving force).






