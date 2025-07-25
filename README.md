# GLoBES-EFT

The phenomenological results presented in this paper were obtained using a custom-built probability engine for the GLoBES framework, which we call GLoBES-EFT. This plugin implements a complete and consistent framework for both Weak Effective Field Theory (WEFT) and Standard Model Effective Field Theory (SMEFT). It correctly models new physics effects in neutrino production, propagation, and detection, and handles the renormalisation group evolution (RGE) using the results in Ref.~\cite{Gonzalez-Alonso:2017iyc} and the matching between the two EFTs given in \cref{sec:matching}. This section serves as a technical manual for users who wish to employ GLoBES-EFT in their own analyses. 

Both production and detection processes are implemented through the use of production coefficients $p_{XY, \alpha}^{S, jk}$ and detection coefficients $d_{XY, \beta}^{\mathrm{Type}, \mathrm{rs}}$, as introduced in \cref{sec:Productions} and \cref{sec:Detection}, respectively. The production coefficients depend on the parent meson $S$, while the detection coefficients depend on the type of cross section. While GLoBES already separates different cross section types, it typically treats fluxes as total fluxes without distinguishing their parent meson origin. However, in our EFT formalism, this separation is essential: each total flux must be decomposed into its contributions from individual parent mesons. To accommodate this, both production and detection coefficients are provided in external files. Each file begins with a first column that specifies the neutrino energy---in units of GeV for production, and in $\log(E)$ for detection. The remaining columns contain the EFT coefficient values, arranged in a fixed and structured order:

- The outermost loop is over the Lorentz structure pairs $XY \in \{\text{LL}, \text{LR}, \text{LS}, \text{LP}, \text{LT}, \dots\}$.
- For each fixed Lorentz structure, the coefficients vary over lepton flavor indices $\alpha = e, \mu, \tau$, first for neutrinos and then for antineutrinos.

To link these files to the fluxes and cross sections in a GLoBES experiment definition, one specifies the production and detection coefficients as follows:

```
nuflux(#DplusNuMode)<
   @flux_file="./FluxesDplusNuModeND.txt"
   @time = 3.5 /* years */
   @power = 11.0  
   @eft_coeff_file = "./plus_ProdALLDpDecay.txt"
   @norm = 4.415e8
   @quark_flavors = 1 : 0
>
```

```
cross(#CC)<
  @cross_file = "./xsec_cc_Oxygen_LQCD_LL_divided_by_16.txt"
  @eft_coeff_file = "./DetAllCCQEOxygenLQCDS1Correct2.txt"
  @quark_flavors = 0 : 0
>
```

The `@eft_coeff_file` field points to the relevant coefficient file, and `@quark_flavors` specifies the quark flavor indices. Once these are properly defined, the GLoBES-EFT engine can be initialized via:

```
smeft_init_probability_engine_3();
```

and you must then register the engine with GLoBES using:

```
glbRegisterProbabilityEngine(MAX_PARAMS, &smeft_probability_matrix,
                                 &smeft_set_oscillation_parameters,
                                 &smeft_get_oscillation_parameters, NULL);
```

where `MAX_PARAMS` is the number of parameters registered within Globes which for our case of having both WEFT and SMEFT corresponds to 2405. As a first check, to verify that coefficients have been successfully loaded, the following functions can be used:

```
double glbEFTFluxCoeff(int exp, int flux_ident, int X, int Y, int alpha, double E)
```

```
double glbEFTXSecCoeff(int exp, int xsec_ident, int X, int Y, int alpha,int cp, double E)
```

these return the values of $p_{XY, \alpha}(E)$ and $d_{XY, \alpha}(E)$, respectively. The `flux_ident` and `xsec_ident` specify the flux and cross section by their order of definition in the GLoBES files, while `cp` is $+1$ for neutrinos and $-1$ for antineutrinos. 

The plugin supports both WEFT and SMEFT modes. If nothing is specified, the code runs in WEFT mode by default. To switch to SMEFT mode, use:

```
 smeft_glbSetOscParamByName(test_values, 1, "SMEFT_FLAG");
```

Note that in SMEFT mode, the Wilson coefficients are assumed to be defined at 1 TeV. The code then runs the RGE evolution to the electroweak scale, performs matching to WEFT, and continues the evolution within WEFT to 2 GeV or to the $b$-quark mass, depending on the operator. In our DUNE implementation, operators involving $b$-quarks are not probed, so all WEFT coefficients are evaluated at 2 GeV. In contrast, WEFT mode assumes the user-defined Wilson coefficients are already at 2 GeV.

For convenience, it is recommended to use `smeft_glbSetOscParamByName` to assign human-readable names to the SMEFT and WEFT parameters as defined and shown in example `example-smeft`. A list of appropriate names is provided in the global array `char smeft_param_strings[][64]`. The naming conventions are detailed in \cref{tab:weft_params_globes} for WEFT and \cref{tab:smeft_params_globes} for SMEFT.

## WEFT Parameter Names

| WEFT Operator | GLoBES Name (Modulus) | GLoBES Name (Argument) |
|---|---|---|
| **Charged Current (CC) Operators** | | |
| $[\epsilon_L^{jk}]_{\alpha\beta}$ | `ABS_EPS_CC_L_jk_pr` | `ARG_EPS_CC_L_jk_pr` |
| $[\epsilon_R^{jk}]_{\alpha\beta}$ | `ABS_EPS_CC_R_jk_pr` | `ARG_EPS_CC_R_jk_pr` |
| $[\epsilon_S^{jk}]_{\alpha\beta}$ | `ABS_EPS_CC_S_jk_pr` | `ARG_EPS_CC_S_jk_pr` |
| $[\epsilon_P^{jk}]_{\alpha\beta}$ | `ABS_EPS_CC_P_jk_pr` | `ARG_EPS_CC_P_jk_pr` |
| $[\epsilon_T^{jk}]_{\alpha\beta}$ | `ABS_EPS_CC_T_jk_pr` | `ARG_EPS_CC_T_jk_pr` |
| **Neutral Current (NC) Operators** | | |
| $[\epsilon_{L}^{m,f}]_{\alpha\beta}$ | `ABS_EPS_NC_L_ff_pr` | `ARG_EPS_NC_L_ff_pr` |
| $[\epsilon_{R}^{m,f}]_{\alpha\beta}$ | `ABS_EPS_NC_R_ff_pr` | `ARG_EPS_NC_R_ff_pr` |
| $[\epsilon_{L}^{m,f}]_{\alpha\alpha}$ | `EPS_NC_L_ff_pp` | --- |
| $[\epsilon_{R}^{m,f}]_{\alpha\alpha}$ | `EPS_NC_R_ff_pp` | --- |

WEFT parameter names available in GLoBES-EFT. Lepton flavor indices $\alpha, \beta$ correspond to GLoBES strings `E, MU, TAU, S1, ...`. Quark indices $j,k$ correspond to `u, c` and `d, s, b` respectively, forming pairs like `ud` or `us`. Matter fermion fields $f$ correspond to `E, u, d`.

## SMEFT Parameter Names

| SMEFT Operator | GLoBES Name (Modulus) | GLoBES Name (Argument) |
|---|---|---|
| $[c_{HWB}]$ | `ABS_C_HWB` | `ARG_C_HWB` |
| $[c_{HD}]$ | `ABS_C_HD` | `ARG_C_HD` |
| $[c_{H l}^{(3)}]_{pr}$ | `ABS_C3_PHI_L_pr` | `ARG_C3_PHI_L_pr` |
| $[c_{H l}^{(1)}]_{pr}$ | `ABS_C_PHI_L_pr` | `ARG_C_PHI_L_pr` |
| $[c_{H e}]_{pr}$ | `ABS_C_PHI_E_pr` | `ARG_C_PHI_E_pr` |
| $[c_{H q}^{(3)}]_{j_q k_q}$ | `ABS_C3_PHI_Q_jq_kq` | `ARG_C3_PHI_Q_jq_kq` |
| $[c_{H q}^{(1)}]_{j_q k_q}$ | `ABS_C_PHI_Q_jq_kq` | `ARG_C_PHI_Q_jq_kq` |
| $[c_{H u}]_{j_u k_u}$ | `ABS_C_PHI_U_ju_ku` | `ARG_C_PHI_U_ju_ku` |
| $[c_{H d}]_{j_d k_d}$ | `ABS_C_PHI_D_jd_kd` | `ARG_C_PHI_D_jd_kd` |
| $[c_{H ud}]_{j_u k_d}$ | `ABS_C_PHI_UD_ju_kd` | `ARG_C_PHI_UD_ju_kd` |
| $[c_{ll}]_{prst}$ | `ABS_C_LL_pr_st` | `ARG_C_LL_pr_st` |
| $[c_{le}]_{prst}$ | `ABS_C_LE_pr_st` | `ARG_C_LE_pr_st` |
| $[c_{lq}^{(3)}]_{pr, j_q k_q}$ | `ABS_C3_LQ_pr_jq_kq` | `ARG_C3_LQ_pr_jq_kq` |
| $[c_{lq}^{(1)}]_{pr, j_q k_q}$ | `ABS_C_LQ_pr_jq_kq` | `ARG_C_LQ_pr_jq_kq` |
| $[c_{lu}]_{pr, j_u k_u}$ | `ABS_C_LU_pr_ju_ku` | `ARG_C_LU_pr_ju_ku` |
| $[c_{ld}]_{pr, j_d k_d}$ | `ABS_C_LD_pr_jd_kd` | `ARG_C_LD_pr_jd_kd` |
| $[c_{ledq}]_{pr, j_d k_q}$ | `ABS_C_LEDQ_pr_jd_kq` | `ARG_C_LEDQ_pr_jd_kq` |
| $[c_{lequ}^{(1)}]_{pr, j_q k_u}$ | `ABS_C_LEQU_pr_jq_ku` | `ARG_C_LEQU_pr_jq_ku` |
| $[c_{lequ}^{(3)}]_{pr, j_q k_u}$ | `ABS_C3_LEQU_pr_jq_ku` | `ARG_C3_LEQU_pr_jq_ku` |

SMEFT operator parameter names available in GLoBES-EFT. Lepton indices $p,r,s,t$ correspond to `E, MU, TAU`. Up-type quark indices $j_u, k_u$ correspond to `u, c, t`; down-type $j_d, k_d$ to `d, s, b`; and quark doublet families $j_q, k_q$ to `Q1, Q2, Q3`.

The user might also want to have access to the oscillation probabilities used. In GLoBES-EFT, these are in fact the pseudo-probability $\tilde{P}_{\alpha \beta}^S\left(E_\nu, L\right)$ defined in \cref{eq:tildeP}. The main difference is that this pseudo-probability depends on both the parent meson and the type of cross section. Instead of using the native GLoBES functions to extract probabilities, one should use the following modified ones:

```
double smeft_glbVacuumProbability(int initial_flavour, int final_flavour,
            int cp_sign, double E, double L,int flux_id, int cross_id);
double smeft_glbConstantDensityProbability(int initial_flavour, int final_flavour,
            int cp_sign, double E, double L, double rho,int flux_id, int cross_id);
double smeft_glbProfileProbability(int exp,int initial_flavour, int final_flavour,
            int panti, double energy,int flux_id, int cross_id);
double smeft_glbFilteredConstantDensityProbability(int exp,int initial_flavour,
            int final_flavour, int panti, double energy,int flux_id, int cross_id);
```

These functions behave like their native counterparts, but require two extra arguments `flux_id` and `cross_id` which specify the flux and cross section by order of definition in the GLB files. 

After using the GLoBES-EFT engine, it is good practice to release the small amount of memory it allocates by calling:

```
smeft_free_probability_engine();
```
