# Immune Responses May Make HIV-1 Therapeutic Interfering Particles Less Effective
This repository contains the code used to simulate the models and generate figures for the manuscript "Immune Responses May Make HIV-1 Therapeutic Interfering Particles Less Effective" (Dodd &amp; de Boer). Note that *grind.R* should be in the same working directory as the scripts in order to simulate the model.

## Scripts:
- *TIP_dynamics_without_immune_response.R*: Simulates the basic TIP model without an immune response (Eqs. 7-14 in the main text) to generate Fig. 2, Fig. 3, and Fig. 4.
- *TIP_dynamics_eclipse_phase_immune_response.R*: Simulates the TIP model with an immune response where effector cells kill virally-infected cells during their eclipse phase. Killing rates of virally-infected cells are either all equal (Eqs. 17-28 in the main text, Fig. 5) or can be independent for WT, mosaic WT, and TIP-infected cells (Eqs. B33-B49 in appendix B.4, Fig. S1)
- *TIP_dynamics_VL_dependent_production_rate.R*: Simulates the basic TIP model without an immune response, except that target cell production becomes a saturation function of the infectous viral load, generating Fig. S2.
- *TIP_dynamics_multiphase_killing.R*: Simulates the extended immune response model where virally-infected cells can be killed during and after eclipse phase (Eqs. B25-B31 in appendix B.3, generating Fig. S4)
- *TIP_sensitivity_analysis.R*: Perform global sensitivity analysis (Johnson's relative weight) of the parameters of the basic TIP model (Eqs. 7-14). Generates Fig. S3.
