# GLoBES-EFT

The phenomenological results presented in this paper were obtained using a custom-built probability engine for the GLoBES framework, which we call \texttt{GLoBES-EFT}. This plugin implements a complete and consistent framework for both Weak Effective Field Theory (WEFT) and Standard Model Effective Field Theory (SMEFT). It correctly models new physics effects in neutrino production, propagation, and detection, and handles the renormalisation group evolution (RGE) using the results in Ref.~\cite{Gonzalez-Alonso:2017iyc} and the matching between the two EFTs given in \cref{sec:matching}. This section serves as a technical manual for users who wish to employ \texttt{GLoBES-EFT} in their own analyses. 

Both production and detection processes are implemented through the use of production coefficients $p_{XY, \alpha}^{S, jk}$ and detection coefficients $d_{XY, \beta}^{\mathrm{Type}, \mathrm{rs}}$, as introduced in \cref{sec:Productions} and \cref{sec:Detection}, respectively. The production coefficients depend on the parent meson $S$, while the detection coefficients depend on the type of cross section. While GLoBES already separates different cross section types, it typically treats fluxes as total fluxes without distinguishing their parent meson origin. However, in our EFT formalism, this separation is essential: each total flux must be decomposed into its contributions from individual parent mesons. To accommodate this, both production and detection coefficients are provided in external files. Each file begins with a first column that specifies the neutrino energy---in units of GeV for production, and in $\log(E)$ for detection. The remaining columns contain the EFT coefficient values, arranged in a fixed and structured order:

\begin{itemize}
    \item The outermost loop is over the Lorentz structure pairs $XY \in \{\text{LL}, \text{LR}, \text{LS}, \text{LP}, \text{LT}, \dots\}$.
    \item For each fixed Lorentz structure, the coefficients vary over lepton flavor indices $\alpha = e, \mu, \tau$, first for neutrinos and then for antineutrinos.
\end{itemize}

To link these files to the fluxes and cross sections in a GLoBES experiment definition, one specifies the production and detection coefficients as follows:

\begin{verbatim}
nuflux(#DplusNuMode)<
   @flux_file="./FluxesDplusNuModeND.txt"
   @time = 3.5 /* years */
   @power = 11.0  
   @eft_coeff_file = "./plus_ProdALLDpDecay.txt"
   @norm = 4.415e8
   @quark_flavors = 1 : 0
>
\end{verbatim}

\begin{verbatim}
cross(#CC)<
  @cross_file = "./xsec_cc_Oxygen_LQCD_LL_divided_by_16.txt"
  @eft_coeff_file = "./DetAllCCQEOxygenLQCDS1Correct2.txt"
  @quark_flavors = 0 : 0
>
\end{verbatim}
The \texttt{@eft\_coeff\_file} field points to the relevant coefficient file, and \texttt{@quark\_flavors} specifies the quark flavor indices. Once these are properly defined, the \texttt{GLoBES-EFT} engine can be initialized via:
\begin{verbatim}
smeft_init_probability_engine_3();
\end{verbatim}

and you must then register the engine with GLoBES using:

\begin{verbatim}
glbRegisterProbabilityEngine(MAX_PARAMS, &smeft_probability_matrix,
                                 &smeft_set_oscillation_parameters,
                                 &smeft_get_oscillation_parameters, NULL);

\end{verbatim}
where \texttt{MAX\_PARAMS} is the number of parameters registered within Globes which for our case of having both WEFT and SMEFT corresponds to 2405. As a first check, to verify that coefficients have been successfully loaded, the following functions can be used:

\begin{verbatim}
double glbEFTFluxCoeff(int exp, int flux_ident, int X, int Y, int alpha, double E)
\end{verbatim}

\begin{verbatim}
double glbEFTXSecCoeff(int exp, int xsec_ident, int X, int Y, int alpha,int cp, double E)
\end{verbatim}
these return the values of $p_{XY, \alpha}(E)$ and $d_{XY, \alpha}(E)$, respectively. The \texttt{flux\_ident} and \texttt{xsec\_ident} specify the flux and cross section by their order of definition in the GLoBES files, while \texttt{cp} is $+1$ for neutrinos and $-1$ for antineutrinos. 

The plugin supports both WEFT and SMEFT modes. If nothing is specified, the code runs in WEFT mode by default. To switch to SMEFT mode, use:

\begin{verbatim}
 smeft_glbSetOscParamByName(test_values, 1, "SMEFT_FLAG");
\end{verbatim}

Note that in SMEFT mode, the Wilson coefficients are assumed to be defined at 1 TeV. The code then runs the RGE evolution to the electroweak scale, performs matching to WEFT, and continues the evolution within WEFT to 2 GeV or to the $b$-quark mass, depending on the operator. In our DUNE implementation, operators involving $b$-quarks are not probed, so all WEFT coefficients are evaluated at 2 GeV. In contrast, WEFT mode assumes the user-defined Wilson coefficients are already at 2 GeV.

For convenience, it is recommended to use \texttt{smeft\_glbSetOscParamByName} to assign human-readable names to the SMEFT and WEFT parameters as defined and shown in example \texttt{example-smeft}. A list of appropriate names is provided in the global array \texttt{char smeft\_param\_strings[][64]}. The naming conventions are detailed in \cref{tab:weft_params_globes} for WEFT and \cref{tab:smeft_params_globes} for SMEFT.


\begin{table}
    \centering
    \begin{tabular}{lll}
        \toprule
        \textbf{WEFT Operator} & \textbf{GLoBES Name (Modulus)} & \textbf{GLoBES Name (Argument)} \\
        \midrule
        \multicolumn{3}{c}{\textit{Charged Current (CC) Operators}} \\
        $[\epsilon_L^{jk}]_{\alpha\beta}$ & \texttt{ABS\_EPS\_CC\_L\_$jk$\_$p r$} & \texttt{ARG\_EPS\_CC\_L\_$jk$\_$p r$} \\
        $[\epsilon_R^{jk}]_{\alpha\beta}$ & \texttt{ABS\_EPS\_CC\_R\_$jk$\_$p r$} & \texttt{ARG\_EPS\_CC\_R\_$jk$\_$p r$} \\
        $[\epsilon_S^{jk}]_{\alpha\beta}$ & \texttt{ABS\_EPS\_CC\_S\_$jk$\_$p r$} & \texttt{ARG\_EPS\_CC\_S\_$jk$\_$p r$} \\
        $[\epsilon_P^{jk}]_{\alpha\beta}$ & \texttt{ABS\_EPS\_CC\_P\_$jk$\_$p r$} & \texttt{ARG\_EPS\_CC\_P\_$jk$\_$p r$} \\
        $[\epsilon_T^{jk}]_{\alpha\beta}$ & \texttt{ABS\_EPS\_CC\_T\_$jk$\_$p r$} & \texttt{ARG\_EPS\_CC\_T\_$jk$\_$p r$} \\
        \midrule
        \multicolumn{3}{c}{\textit{Neutral Current (NC) Operators}} \\
        $[\epsilon_{L}^{m,f}]_{\alpha\beta}$ & \texttt{ABS\_EPS\_NC\_L\_$f f$\_$p r$} & \texttt{ARG\_EPS\_NC\_L\_$f f$\_$p r$} \\
        $[\epsilon_{R}^{m,f}]_{\alpha\beta}$ & \texttt{ABS\_EPS\_NC\_R\_$f f$\_$p r$} & \texttt{ARG\_EPS\_NC\_R\_$f f$\_$p r$} \\
        $[\epsilon_{L}^{m,f}]_{\alpha\alpha}$ & \texttt{EPS\_NC\_L\_$f f$\_$p p$} & --- \\
        $[\epsilon_{R}^{m,f}]_{\alpha\alpha}$ & \texttt{EPS\_NC\_R\_$f f$\_$p p$} & --- \\
        \bottomrule
    \end{tabular}

    \caption{WEFT parameter names available in \texttt{GLoBES-EFT}. Lepton flavor indices $\alpha, \beta$ correspond to GLoBES strings \texttt{E, MU, TAU, S1, ...}. Quark indices $j,k$ correspond to \texttt{u, c} and \texttt{d, s, b} respectively, forming pairs like \texttt{ud} or \texttt{us}. Matter fermion fields $f$ correspond to \texttt{E, u, d}.}
    \label{tab:weft_params_globes}
\end{table}



\begin{table}
    \centering 
    \begin{tabular}{lll}
        \toprule
        \textbf{SMEFT Operator} & \textbf{GLoBES Name (Modulus)} & \textbf{GLoBES Name (Argument)} \\
        \midrule
        $[c_{HWB}]$ & \texttt{ABS\_C\_HWB} & \texttt{ARG\_C\_HWB} \\
        $[c_{HD}]$ & \texttt{ABS\_C\_HD} & \texttt{ARG\_C\_HD} \\
        $[c_{H l}^{(3)}]_{pr}$ & \texttt{ABS\_C3\_PHI\_L\_$pr$} & \texttt{ARG\_C3\_PHI\_L\_$pr$} \\
        $[c_{H l}^{(1)}]_{pr}$ & \texttt{ABS\_C\_PHI\_L\_$pr$} & \texttt{ARG\_C\_PHI\_L\_$pr$} \\
        $[c_{H e}]_{pr}$ & \texttt{ABS\_C\_PHI\_E\_$pr$} & \texttt{ARG\_C\_PHI\_E\_$pr$} \\
        $[c_{H q}^{(3)}]_{j_q k_q}$ & \texttt{ABS\_C3\_PHI\_Q\_$j_q k_q$} & \texttt{ARG\_C3\_PHI\_Q\_$j_q k_q$} \\
        $[c_{H q}^{(1)}]_{j_q k_q}$ & \texttt{ABS\_C\_PHI\_Q\_$j_q k_q$} & \texttt{ARG\_C\_PHI\_Q\_$j_q k_q$} \\
        $[c_{H u}]_{j_u k_u}$ & \texttt{ABS\_C\_PHI\_U\_$j_u k_u$} & \texttt{ARG\_C\_PHI\_U\_$j_u k_u$} \\
        $[c_{H d}]_{j_d k_d}$ & \texttt{ABS\_C\_PHI\_D\_$j_d k_d$} & \texttt{ARG\_C\_PHI\_D\_$j_d k_d$} \\
        $[c_{H ud}]_{j_u k_d}$ & \texttt{ABS\_C\_PHI\_UD\_$j_u k_d$} & \texttt{ARG\_C\_PHI\_UD\_$j_u k_d$} \\
        $[c_{ll}]_{prst}$ & \texttt{ABS\_C\_LL\_$pr$\_$st$} & \texttt{ARG\_C\_LL\_$pr$\_$st$} \\
        $[c_{le}]_{prst}$ & \texttt{ABS\_C\_LE\_$pr$\_$st$} & \texttt{ARG\_C\_LE\_$pr$\_$st$} \\
        $[c_{lq}^{(3)}]_{pr, j_q k_q}$ & \texttt{ABS\_C3\_LQ\_$pr$\_$j_q k_q$} & \texttt{ARG\_C3\_LQ\_$pr$\_$j_q k_q$} \\
        $[c_{lq}^{(1)}]_{pr, j_q k_q}$ & \texttt{ABS\_C\_LQ\_$pr$\_$j_q k_q$} & \texttt{ARG\_C\_LQ\_$pr$\_$j_q k_q$} \\
        $[c_{lu}]_{pr, j_u k_u}$ & \texttt{ABS\_C\_LU\_$pr$\_$j_u k_u$} & \texttt{ARG\_C\_LU\_$pr$\_$j_u k_u$} \\
        $[c_{ld}]_{pr, j_d k_d}$ & \texttt{ABS\_C\_LD\_$pr$\_$j_d k_d$} & \texttt{ARG\_C\_LD\_$pr$\_$j_d k_d$} \\
        $[c_{ledq}]_{pr, j_d k_q}$ & \texttt{ABS\_C\_LEDQ\_$pr$\_$j_d k_q$} & \texttt{ARG\_C\_LEDQ\_$pr$\_$j_d k_q$} \\
        $[c_{lequ}^{(1)}]_{pr, j_q k_u}$ & \texttt{ABS\_C\_LEQU\_$pr$\_$j_q k_u$} & \texttt{ARG\_C\_LEQU\_$pr$\_$j_q k_u$} \\
        $[c_{lequ}^{(3)}]_{pr, j_q k_u}$ & \texttt{ABS\_C3\_LEQU\_$pr$\_$j_q k_u$} & \texttt{ARG\_C3\_LEQU\_$pr$\_$j_q k_u$} \\
        \bottomrule
    \end{tabular}

    \caption{SMEFT operator parameter names available in \texttt{GLoBES-EFT}. Lepton indices $p,r,s,t$ correspond to \texttt{E, MU, TAU}. Up-type quark indices $j_u, k_u$ correspond to \texttt{u, c, t}; down-type $j_d, k_d$ to \texttt{d, s, b}; and quark doublet families $j_q, k_q$ to \texttt{Q1, Q2, Q3}.}
    \label{tab:smeft_params_globes}
\end{table}
The user might also want to have access to the oscillation probabilities used. In \texttt{GLoBES-EFT}, these are in fact the pseudo-probability $\tilde{P}_{\alpha \beta}^S\left(E_\nu, L\right)$ defined in \cref{eq:tildeP}. The main difference is that this pseudo-probability depends on both the parent meson and the type of cross section. Instead of using the native GLoBES functions to extract probabilities, one should use the following modified ones:

\begin{verbatim}
double smeft_glbVacuumProbability(int initial_flavour, int final_flavour,
            int cp_sign, double E, double L,int flux_id, int cross_id);
double smeft_glbConstantDensityProbability(int initial_flavour, int final_flavour,
            int cp_sign, double E, double L, double rho,int flux_id, int cross_id);
double smeft_glbProfileProbability(int exp,int initial_flavour, int final_flavour,
            int panti, double energy,int flux_id, int cross_id);
double smeft_glbFilteredConstantDensityProbability(int exp,int initial_flavour,
            int final_flavour, int panti, double energy,int flux_id, int cross_id);
\end{verbatim}

These functions behave like their native counterparts, but require two extra arguments \texttt{flux\_id} and \texttt{cross\_id} which specify the flux and cross section by order of definition in the GLB files. 

After using the \texttt{GLoBES-EFT} engine, it is good practice to release the small amount of memory it allocates by calling:

\begin{verbatim}
smeft_free_probability_engine();
\end{verbatim}
