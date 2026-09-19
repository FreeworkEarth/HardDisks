# Understanding HardDisks: pressure, sound, Knudsen number, work and predictive control

Date: 2026-09-10, local session date. Existing project notes also use 260911/260912 prefixes.

Status: explanatory research note, not a completed-results report or an instruction to run experiments.

Prepared at Chris's request. Existing notes and selected code were inspected. No simulations, process termination, production-code changes, commits or independent numerical reanalysis were performed. Numerical results quoted from Chris's pasted Claude report remain attributed to that report. Other agents are updating some files concurrently; this note does not certify their latest outputs.

## Start here: what are the three papers trying to establish?

1. **Paper 1 — measurement and validation:** Does the simulator reproduce equilibrium properties, and when does the divider measure bulk sound speed rather than a finite-box response?
2. **Paper 2 — driven thermodynamics:** Where does piston work go, and how much exceeds the appropriate reversible reference?
3. **Possible Paper 3 — predictive control:** How much information about the recent past must a controller retain to achieve reliable energy transfer or recovery?

These are distinct questions. Success in one supports the next but does not prove it.

> Pressure is the push now. Compressibility is how the push changes when squeezed. Sound speed is how quickly a small compression propagates.

## 1. Geometry: low density does not mean only two particles

Let N_s be the disks in one compartment, r their radius, sigma=2r their diameter, H the height and L_0 the compartment length. With the campaign's nominal geometric area:

$$
A=L_0H,\qquad n=N_s/A,\qquad
\eta=n\pi r^2=\frac{N_s\pi r^2}{L_0H}.
$$

On route A, N_s=50, r=0.5 and H=10:

$$
\eta=\frac{3.92699}{L_0}.
$$

| Packing fraction eta | L_0 in disk diameters | Disks per compartment |
|---|---:|---:|
| 0.02 | 196.35 | 50 |
| 0.10 | 39.27 | 50 |
| 0.65 | 6.04 | 50 |

At eta=0.02 the same 50 disks occupy a longer compartment; there are not just two disks.

Three experiments must not be conflated:

- Fixed particle number and radius, changing area: changes density.
- Fixed density, radius and shape, increasing particle number and area together: tests system size.
- Fixed box, adding particles: changes density and particle number together.

Shrinking disks while adding particles also enlarges the box in units of disk diameter. Comparing that with a fixed-diameter campaign requires matching other dimensionless parameters, not just eta.

Physical walls exclude disk centres from a boundary layer. Nominal area, centre-accessible area and effective length in a divider-mode formula are not automatically interchangeable. State the convention and use it consistently.

## 2. Pressure, Z and actual compressibility

In 2D, wall pressure is normal momentum impulse per wall length per time:

$$
P_{\rm wall}=\frac{\sum_j\Delta p_{\perp,j}}{H\,\Delta t}.
$$

Its units are force/length, unlike the force/area of a 3D fluid. For an ideal gas P=n k_BT. Hard-disk nonideality is summarized by

$$
Z=\frac{P}{nk_BT}=\frac{PA}{Nk_BT}.
$$

Z=1 means ideal-gas pressure at the same n and T. Z=8.4 means 8.4 times that pressure. Despite its name, the “compressibility factor” Z is a pressure ratio.

Actual isothermal compressibility is

$$
\kappa_T=-\frac1A\left(\frac{\partial A}{\partial P}\right)_{T,N}.
$$

Large kappa means easy to squeeze; its inverse is a stiffness. A single pressure value does not tell us how pressure changes on compression.

For hard disks, P=n k_BT Z(eta), with eta proportional to n:

$$
\left(\frac{\partial P}{\partial n}\right)_T
=k_BT[Z+\eta Z'],\qquad
\kappa_T=\frac{1}{nk_BT[Z+\eta Z']}.
$$

An adiabatic compression additionally heats the gas. For the homogeneous equilibrium 2D hard-disk fluid, U=Nk_BT and dU=-P dA imply

$$
\frac{dT}{T}=Z\frac{dn}{n}.
$$

Using mass density rho=mn, the hydrodynamic adiabatic sound speed is

$$
\boxed{
c_s^2=\left(\frac{\partial P}{\partial\rho}\right)_s
=\frac{1}{\rho\kappa_S}
=\frac{k_BT}{m}[Z+\eta Z'+Z^2].
}
$$

Here s is entropy per unit mass. The three contributions are: more particles per area; changing excluded-volume pressure; and compression heating. For an ideal gas Z=1, Z'=0, so c_s=sqrt(2k_BT/m), or sqrt(2) in the normalized units.

This is a bulk-fluid relation. Applying it to a small, narrow, dynamically driven box is a physical approximation to test.

## 3. Why pressure checks support a sound-speed paper

Pressure tests the static equation of state and momentum transfer. Divider oscillations test dynamic response, moving boundaries, the mode formula and frequency inference. Conservation tests numerical bookkeeping.

A correct-looking number can hide compensating errors. The reported old temperature initialization gives mean T_i near 0.98 rather than 1:

$$
\sqrt{0.98}=0.98995.
$$

That lowers a hard-disk frequency by about 1%. Correcting temperature can reveal a previously concealed positive offset. This is a reason to use complementary tests, not to discard all prior data.

Neither pressure nor sound agreement proves all simulator physics. Wall and collision-based pressure also share dynamics and conservation constraints, so they are not automatically statistically independent checks.

## 4. What the frequency experiment measures

The measured quantity is the divider trajectory x(t). A model converts its oscillation frequency to a sound-speed estimate:

$$
\nu=c_s\,q(M,\text{geometry}).
$$

The documented Roman-style relation uses q=K/(2 pi L_eff) and cot(K)=[M/(2N_s m)]K. Its effective-length convention must match the theory and geometry; it is not a universal formula.

Several divider masses test whether one c_s explains several oscillation frequencies.

- **25 repeats:** independent realizations for ensemble variability.
- **25 target oscillations:** intended duration of each realization, not a guarantee of 25 observed clean periods.
- **FFT bin spacing:** approximately 1/T_obs.
- **Sub-bin fitting:** can locate a frequency between FFT bins, subject to model and noise assumptions.

Chris's latest pasted report says damped-cosine fitting reached 0.02% fit precision and a temperature slope 1.0 +/- 0.1 with answer-independent selection. That is useful new evidence. My earlier warning against promising precision in advance did not deny that such an outcome was possible.

I found [fit_nu_damped.py](../../hspist3/validation/fit_nu_damped.py), which obtains sigma_nu from the local parameter covariance returned by curve_fit. I did not independently reproduce the reported aggregate precision, slope or quality-cut statistics.

Distinguish **fit precision** from **total physical accuracy**. Correlated residuals, mode mismatch, geometry and between-run scatter are not automatically included in that covariance. SciPy explicitly describes it as approximate. [SciPy curve_fit documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html)

One targeted uncertainty-calibration and selection-coverage check is sufficient; this is not a request for an endless audit.

## 5. Knudsen number: physical intuition and the corrected calculation

The mean free path lambda is the typical distance travelled between particle-particle collisions. It is not the disk spacing, disk diameter, sound wavelength or distance between wall bounces.

The Knudsen number compares that microscopic distance with a relevant macroscopic length:

$$
\mathrm{Kn}_L=\frac{\lambda}{L}.
$$

Collisions provide opportunities to redistribute momentum and energy locally. Many collisions over the scale on which the macroscopic state varies make a fluid description more plausible. With few collisions, individual flight histories matter.

Small Kn supports, but does not alone guarantee, hydrodynamics. There is no universal particle-count threshold separating a fluid from a kinetic regime.

### 5.1 Where 0.278 comes from

For diameter sigma, a collision occurs for impact parameters between -sigma and +sigma: a strip of width 2sigma. With isotropic Maxwell velocities:

$$
\nu_c=2n\sigma g(\sigma)\langle v_{\rm rel}\rangle,\qquad
\langle v_{\rm rel}\rangle=\sqrt2\langle v\rangle,
$$

$$
\lambda=\frac{\langle v\rangle}{\nu_c}
=\frac{1}{2\sqrt2\,n\sigma g(\sigma)}.
$$

Here g(sigma) describes contact correlations and tends to 1 in the dilute limit. Since n=4eta/(pi sigma^2):

$$
\boxed{\lambda_{\rm dilute}
=\frac{\pi}{8\sqrt2}\frac{\sigma}{\eta}
\simeq0.27768\frac{\sigma}{\eta}.}
$$

The earlier 0.555 coefficient missed the collision-strip factor. This convention agrees with the mean-free-path expression in the cited Enskog study. [Non-equilibrium dynamics of dense gas under tight confinement, section 2](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/nonequilibrium-dynamics-of-dense-gas-under-tight-confinement/3E60877A4D95AFCB96DAFF8AA00F3062)

### 5.2 There is more than one relevant scale

At eta=0.02, the dilute estimates give

$$
\lambda\simeq13.88\sigma,\quad L_0\simeq196.35\sigma,\quad
\mathrm{Kn}_{L_0}\simeq0.0707,\quad
\mathrm{Kn}_{H}\simeq1.39.
$$

The same compartment is long compared with a free path in one direction but narrow in the other. Specular wall bounces do not thermalize velocities as a thermal wall would.

For a sound mode also consider wave number k and angular frequency omega:

$$
k\lambda,\qquad \omega\tau_c,\qquad
\tau_c=1/\nu_c=\lambda/\langle v\rangle.
$$

These compare flight distance with spatial variation and collision time with temporal variation. Do not confuse nu in cycles/time with omega=2 pi nu.

Along route A, lambda/L_0 is approximately constant only in the dilute approximation; contact correlations change at finite density. Because H stays fixed, lambda/H is not constant.

### 5.3 Does halving lambda prove a fourfold-smaller sound-speed offset?

Only conditionally. If a justified model has

$$
\frac{c_{\rm eff}-c_s}{c_s}=C(\omega\tau_c)^2+\cdots
$$

and the only erroneous input was tau_c, halving tau_c quarters this term. But C, its sign, regime and boundary dependence still need justification. “A few tenths of a percent” is a scale estimate, not a measured decomposition of the offset.

famB is useful but not a binary verdict. Wall corrections, mode-model errors and other finite-size effects can also shrink with N. A shrinking offset establishes size dependence, not uniquely a kinetic cause.

Compare physically motivated scaling and track k lambda, omega tau_c, aspect ratio and M/(N_s m). Using the same absolute divider masses across sizes changes that last ratio.

### 5.4 The sqrt(3) limit: the legitimate argument and its boundary

There is a real mechanical argument. Under slow uniaxial compression of an ideal collisionless gas with specular boundaries, the single-particle adiabatic invariant gives |p_x|L approximately constant. Thus T_x is proportional to L^-2, and P_x to L^-3 at fixed height and particle number:

$$
\left(\frac{\partial P_x}{\partial\rho}\right)_{\rm directional}
=\frac{3k_BT_x}{m}.
$$

This motivates a directional stiffness speed sqrt(3k_BT_x/m). It is not ordinary locally equilibrated sound, nor a guarantee that a finite-frequency divider mode approaches precisely that number.

My earlier caveat concerned that identification and the “ten particles” threshold, not the existence of the limiting argument.

## 6. Finite-size fits and the proposed next CC prompt

At fixed density and shape, area grows as N and perimeter as sqrt(N). A finite boundary-layer thickness motivates a leading 1/sqrt(N) correction. A 1/N form is a useful sensitivity comparison, not automatically equally well motivated.

$$
Z(N)=Z_\infty+aN^{-1/2},\qquad
\chi^2=\sum_i\frac{[Z_i-Z_{\rm fit}(N_i)]^2}{\sigma_i^2}.
$$

Four independent size estimates and two fitted coefficients give two degrees of freedom. Correlated estimates require a covariance-matrix treatment.

The pasted CC prompt is broadly sensible, with these bounded refinements:

1. **Seed pad:** a diameter-scaled inward placement margin is reasonable, but 1e-3 d is not a universally proven value. Check stored-coordinate wall clearance and pair non-overlap at the largest relevant coordinate scale. Do not change physical walls, radii, particle count, area or eta. Preserve prior artifacts as superseded rather than overwrite them.
2. **Regression gate:** “new nu agrees within its fit error” uses the wrong uncertainty for chaotic finite-duration trajectories. Compare ensemble means using paired seed differences when appropriate, with a predeclared confidence/equivalence criterion including run-to-run variability. Numerical invariant checks remain separate. Ten clean trajectories test that cell, not every geometry.
3. **Pressure forms:** report both fits, residuals and acceptability. Their estimates share data, so errors are correlated. A spread between defensible models is a sensitivity estimate, not automatically a calibrated Gaussian systematic error. Do not give an unacceptable alternative fit equal weight merely to enlarge uncertainty.
4. **Slow work:** the mean at u<=0.02 is a low-speed average, not automatically W_0. A linear extrapolation tests a different assumption. Call the mean a plateau estimate only if residual speed dependence is unresolved at the required accuracy.
5. **GUI:** a separate convenience change, not physical validation. Preserve CLI semantics and invalid-geometry checks; it need not block the scientific report.

Avoid a robust “5.5 sigma below KR” statement without reference and systematic uncertainties. At eta=0.67/0.69 say “these data do not support a reliable bulk extrapolation with this model,” not “no bulk value exists.”

I did not execute the pasted CC prompt.

## 7. Work and the reversible reference

For work positive into the system:

$$
\Delta E_{\rm gas}+\Delta E_{\rm walls}+\Delta E_{\rm spring}
=W_{\rm external}+Q_{\rm in}.
$$

Elastic dynamics do not destroy energy. “Dissipation” requires specifying what organized energy becomes unavailable under which conditions.

For reversible compression:

$$
W_{\rm rev}=-\int_{A_i}^{A_f}P(A,T(A))\,dA.
$$

For a 2D ideal gas initially at T_i:

$$
W_{\rm rev}^{\rm ad}=Nk_BT_i\left(\frac{A_i}{A_f}-1\right),\qquad
W_{\rm rev}^{\rm iso}=Nk_BT_i\ln\frac{A_i}{A_f}.
$$

At N=50, k_BT_i=1 and A_f=0.9 A_i these are 5.56 and 5.27. Their difference is not dissipation: they describe different reversible paths.

Rescaling an entire work reference by an initial-pressure ratio assumes a corresponding correction along the path. Even a constant multiplicative correction to Z changes the adiabatic temperature path through d ln T=Z d ln n. A single hold-time pressure ratio is a heuristic estimate, not a derivation of the finite-box adiabat.

For D=W_0-W_qs:

$$
\operatorname{Var}(D)=\operatorname{Var}(W_0)+
\operatorname{Var}(W_{\rm qs})-2\operatorname{Cov}(W_0,W_{\rm qs}).
$$

Shared temperatures or seeds can create covariance. Uncalibrated reference-model uncertainty belongs beside, not hidden inside, a statistical error bar.

## 8. Jarzynski: an ensemble statement about work

Define beta_th=1/(k_BT), with initial canonical equilibrium and a prescribed protocol:

$$
\boxed{\left\langle e^{-\beta_{\rm th}W}\right\rangle
=e^{-\beta_{\rm th}\Delta F_{\rm eq}}.}
$$

Equivalently,

$$
\Delta F_{\rm eq}=-k_BT\ln\left\langle e^{-\beta_{\rm th}W}\right\rangle,
\qquad \langle W\rangle\geq\Delta F_{\rm eq}.
$$

Brackets average repeated realizations: this is not the exponential of mean work. Rare low-work trajectories can dominate, making finite-sample estimates difficult. The initial ensemble must be canonical; the system need not stay equilibrated during driving. Hamiltonian driving after canonical preparation is possible.

For hard disks, prepare equilibrium allowed positions and the appropriate thermal momentum distribution, including any explicitly chosen constraints. Rescaling each trajectory to exactly the same kinetic energy is not ordinary canonical sampling.

This relation estimates a free-energy difference at T, not the isolated adiabatic reference above. Feedback needs a modified treatment. [Jarzynski, A nonequilibrium equality for free energy differences](https://arxiv.org/abs/cond-mat/9610209)

## 9. Still: memory that predicts versus memory that does not

In Still et al.'s externally driven, no-feedback, heat-bath/Markov setup, X_t is the signal and S_t the responding state. With natural logarithms:

$$
I(U;V)=\sum_{u,v}p(u,v)\ln\frac{p(u,v)}{p(u)p(v)},
$$

$$
I_{\rm mem}(t)=I(S_t;X_t),\qquad
I_{\rm pred}(t)=I(S_t;X_{t+1}).
$$

Their work-step identity and protocol bound are

$$
\beta_{\rm th}\langle W_{\rm diss}^{\,\text{work step }t}\rangle
=I_{\rm mem}(t)-I_{\rm pred}(t),
$$

$$
\sum_t[I_{\rm mem}(t)-I_{\rm pred}(t)]
\leq\beta_{\rm th}\langle W_{\rm diss}^{\rm total}\rangle
\leq\beta_{\rm th}\langle W_{\rm ex}^{\rm total}\rangle.
$$

The definitions distinguish

$$
F_{\rm neq}=F_{\rm eq}+k_BT D_{\rm KL}(p\Vert p_{\rm eq}),\quad
\langle W_{\rm diss}\rangle=\langle W\rangle-\Delta F_{\rm neq},\quad
\langle W_{\rm ex}\rangle=\langle W\rangle-\Delta F_{\rm eq}.
$$

Useful retained information predicts the next signal. Final nonequilibrium free energy prevents identifying all excess work with already dissipated work. In bits, multiply information terms by ln(2). These equations do not automatically apply to an isolated feedback-controlled piston. [Still et al., Thermodynamics of Prediction, equations 9–18](https://arxiv.org/html/1203.3271v2)

## 10. How to map the information question onto this project

The following is my proposed experimental mapping, not a measured result:

- **External signal:** a stochastic piston command with known statistics and correlation time.
- **Responding state:** gas/divider observables after a specified response interval.
- **Prediction target:** a future command or a physically meaningful future outcome.
- **System boundary:** declare gas, bath, actuator and controller.
- **Timing:** distinguish changing the command from the subsequent response.

A deterministic waveform repeated identically has no trial-to-trial input entropy at fixed time; it is not by itself an information-transmission experiment. Pooling different times can merely encode the clock. Randomize messages/protocols explicitly.

A reduced observable is not automatically a sufficient thermodynamic state. Coarse-graining two mutual informations separately does not automatically preserve their difference as a lower bound. Nor is any arbitrary chosen difference guaranteed nonnegative. Check the process assumptions instead of clipping estimates to zero.

Train encoders and select hyperparameters on training runs; evaluate on independent held-out runs and account for temporal dependence. A predictor seeing future samples is not a causal controller.

## 11. The IB package: what exists and what remains

Inspected:

- [Package README](../../IB_Package/IB_DIB_PACKAGE/README.md)
- [Core algorithm](../../IB_Package/IB_DIB_PACKAGE/src/ib/core/algorithm.py)
- [Core configuration](../../IB_Package/IB_DIB_PACKAGE/src/ib/core/config.py)
- [Online module placeholder](../../IB_Package/IB_DIB_PACKAGE/src/ib/online/__init__.py)

The core accepts a discrete joint distribution p(x,y), shape (n_y,n_x), and constructs a compressed variable M:

$$
\max\left[I(Y;M)-\frac1\tau I(X;M)\right].
$$

It implements deterministic annealing with split/merge machinery. The inspected online module is empty, and the README places dynamical IB in a future version. Tests are documented there; I did not run them for this explanatory task.

Deterministic annealing, deterministic IB and dynamical IB are different concepts. Directory naming does not demonstrate all three implementations.

The current core can support an **offline predictive-IB pilot** by constructing X from past observation windows and Y from future outcomes. That still requires data processing and validation. A recurrent online memory and closed-loop controller do not follow automatically.

The package exports information traces in bits; the thermal formulas above use nats. Its annealing tau is not physical time, collision time or automatically a bath-temperature ratio.

## 12. A focused, potentially interesting Paper 3

My assessment: this is scientifically worthwhile if the contribution is a tested physical tradeoff, not just “we applied IB to pistons.”

> At fixed signal reliability or control performance, how much predictive memory is needed, and what energy cost does that choice incur?

A possible causal architecture:

~~~text
past sensor observations + previous memory
                  |
                  v
          compressed memory M_t
                  |
                  v
         bounded piston action a_t
                  |
                  v
     gas/wave/load dynamics and energy ledger
                  |
                  v
           next observation
~~~

For example, encode recent divider position, velocity and pressure into a small number of memory states. Use those states to decide whether a receiving piston should absorb, reflect or return an incoming disturbance.

A predictive representation objective is

$$
\min_{p(m|h)}
I(H_t;M_t)-\beta_{\rm IB}I(M_t;Y_{t+\Delta}),
$$

where H_t is observation history and Y a chosen future target. This is a design objective, not a thermodynamic identity. beta_IB is not beta_th. [Original information bottleneck method](https://arxiv.org/abs/physics/0004057)

Prediction alone may preserve irrelevant detail. A task-aligned control formulation could be

$$
\min_{q,\pi}\;
\mathbb E[C_{\rm task}+a\,W_{\rm cost}]
+b\,I(H_t;M_t),
$$

or a constrained problem fixing task quality and memory budget. Define units and the work-cost functional. Mutual information is not automatically the controller's hardware energy consumption.

A recurring memory could use q(M_t|M_{t-1},O_t), followed by policy pi(a_t|M_t). How to charge memory capacity versus newly acquired information is a modeling decision.

### A small first experiment

1. Fix one validated density, geometry and manageable particle count.
2. Drive one boundary with a randomized, temporally correlated signal. Measure the receiver without adaptive feedback first.
3. Choose one task: decode a message, or capture specified energy in a load by a deadline.
4. Compare a passive baseline, optimized open-loop protocol, simple phase/velocity feedback, and compressed-memory control. A full-observation controller is another benchmark, not automatically a mathematical optimum.
5. Match actuator amplitude, speed, bandwidth, duration and training budgets.
6. Sweep memory budget and input correlation time; test on new seeds and unseen protocols.

Use two controlled pistons after the one-controller question is interpretable. Reflection already returns a wave; a second actuator is useful if coordinated actions provide a distinct, tested advantage.

### What would make the results informative?

- A reproducible tradeoff among information budget, task error and work cost.
- Identification of retained variables: phase, velocity, arrival time or pressure imbalance.
- A change in useful memory when input correlation time crosses propagation or relaxation time.
- A predictive controller outperforming equally constrained nonpredictive controls.
- A negative result showing a simple physical state estimator is as good as IB.

Measure signal fidelity and energy transfer separately. A clear signal need not carry much energy; substantial energy can arrive without preserving a message. For repeated exchanges, distinguish successive messages from reflections and common-driver correlations.

Record work from **both** actuators, load energy, gas energy and initial/final stored energy. A successful receiver may be injecting its own work. Repeated operation also needs a heating/cooling account; a closed elastic gas is not an unlimited equilibrium reservoir.

## 13. Feedback changes the thermodynamic question

Once measurements determine future piston actions, the drive is no longer independent of system response.

In a suitable single-measurement feedback setup, an illustrative generalized relation is

$$
\left\langle e^{-\beta_{\rm th}(W-\Delta F)-i}\right\rangle=1,\qquad
\langle W\rangle\geq\langle\Delta F\rangle-k_BT I,
$$

where i=ln[p(y|s)/p(y)] is trajectory-level measurement information and I=E[i]. This is not a plug-in formula for arbitrary repeated feedback: protocol, support and ensemble assumptions matter. [Sagawa and Ueda, Generalized Jarzynski Equality under Nonequilibrium Feedback Control](https://arxiv.org/abs/0907.4914)

Classical Hamiltonian feedback treatments also exist. [Sagawa, Hamiltonian Derivations](https://arxiv.org/abs/1105.5888)

The practical choice is not “put a thermostat into everything.” Choose a formulation matching the experiment, then define the ledger and ensemble consistently. These are future design decisions, not changes authorized by this note.

## 14. Axons and slime mold: useful analogies, different mechanisms

### Axons

Axons propagate regenerative electrochemical signals involving membrane capacitance, voltage-dependent ionic conductances and ion gradients. A passive hard-disk pressure wave lacks those mechanisms. [Hodgkin and Huxley, 1952](https://physoc.onlinelibrary.wiley.com/doi/10.1113/jphysiol.1952.sp004764)

The useful comparison is functional: delay, noise, bandwidth, energy cost and information needed for a response. “Minimal physical communication channel” is defensible. “Axon simulation” requires a biophysical mapping and additional dynamics.

### Physarum / slime mold

Experimental work links Physarum signal propagation to flow-transported chemical signaling coupled to tube contractions. Its adaptive network is not simply an acoustic channel. [Alim et al., Mechanism of signal propagation in Physarum polycephalum](https://doi.org/10.1073/pnas.1618114114)

A later analogy could examine how geometry and active boundaries influence delayed noisy communication. A biological model would additionally need viscous flow, signal advection/diffusion, active contraction and network remodeling.

Start with the general physical principle. Transfer to biology only when a dimensionless comparison or reduced mathematical model establishes what is shared. Optimality under our chosen objective does not establish that evolution or slime mold optimizes that same objective.

## 15. Is it new? What can we honestly claim?

IB, prediction, energetic efficiency and feedback control already have substantial literature. Susanne Still's [Information theoretic approach to interactive learning](https://arxiv.org/abs/0709.1948) addresses action-dependent learning. [Path Integral Bottleneck work](https://arxiv.org/abs/2505.09896) also connects IB and control. These are starting points, not an exhaustive novelty search.

Potentially distinctive work here combines a mechanically explicit many-particle channel, controlled propagation/confinement, measured work accounting, and a causal memory-limited controller with strong baselines.

Publishable novelty depends on the result and a focused literature comparison. It is not guaranteed by combining familiar topics. Discovering conditions where predictive compression does **not** improve physical efficiency could also be informative.

## 16. What the sources are for — must we repeat their experiments?

“Primary source” means an original research article, not a mandatory experiment.

| Source | Purpose here | Must we repeat its experiment? |
|---|---|---|
| Enskog confinement study, section 5 | Collision-length convention and competing microscopic/channel scales | No. It studies driven channel flow, not our divider experiment. |
| [Kolafa–Rottner equation of state](https://arxiv.org/abs/cond-mat/0608356) | Simulation-based bulk pressure benchmark | No full reproduction. Compare compatible quantities, conventions and uncertainties. |
| Jarzynski, section 8 | Work/free-energy ensemble relation | No original apparatus required; our ensemble and protocol must meet its assumptions. |
| Still et al., section 9 | Precise information/dissipation definitions | No identical apparatus; a claimed test must implement the relevant setup. |
| Original IB, section 12 | Information-preserving compression objective | No. Validate our encoder and generalization. |
| Sagawa–Ueda, section 13 | Feedback changes fluctuation relations | Relevant when actions use measurements. |
| Hodgkin–Huxley and Alim et al., section 14 | Biological mechanisms distinguishing analogy from model | No biological replication required for a toy-model physics paper. |

A bounded optional mean-free-path check could use a dilute equilibrium run, count only particle-particle events, and compute

$$
\nu_{c,\rm measured}=\frac{2C_{\rm pair}}{N\,\Delta t},\qquad
\lambda_{\rm measured}=\frac{\langle v\rangle}{\nu_{c,\rm measured}}.
$$

The 2 counts both participants in each pair collision. Exclude wall events. Match temperature, density convention and finite-boundary assumptions. This is a focused diagnostic, not a request to launch another programme.

For a kinetic explanation of the sound offset, controlled dependence on the relevant dimensionless frequency/length scales is stronger evidence than agreement with a mean-free-path formula alone.

## 17. Reading order and meeting summary

Read sections 1–5 first, then 7–9. Treat 10–15 as a proposed research direction, not completed work.

For Susanne:

> We are separating bulk thermodynamics, finite-box response and numerical error. The corrected mean-free-path estimate supports a kinetic-effects hypothesis but does not quantitatively explain the sound offset yet. We are establishing the correct reversible-work reference before interpreting excess work. A possible next paper would test the memory–performance–energy tradeoff of a causal predictive controller in this mechanically explicit channel, with biological relevance framed as a hypothesis rather than an assumed equivalence.

## Local provenance

Campaign context came from [260911_handoff_state_COWORK.md](260911_handoff_state_COWORK.md), [260911_advisor_meeting_prep_and_commands_COWORK.md](260911_advisor_meeting_prep_and_commands_COWORK.md), previously inspected explanatory notes, and Chris's latest pasted reply. This note supersedes neither accepted trajectories nor analysis outputs.

The reported 0.02% frequency-fit precision and slope were not independently recomputed here. Existing PDF/TeX notes were not edited by this task. Only this new Markdown was created; nothing was staged or committed.
