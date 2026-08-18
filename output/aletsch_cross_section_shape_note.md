# Aletsch cross-section shape: observations, interpolation and power law fit

We first consider the airborne GPR cross-sections acquired during the 2009 and 2011 campaigns. From this dataset, we retain the flow-transverse profiles for which a bedrock reflection can be identified across the full glacier width. Profiles that are strongly oblique or approximately parallel to glacier flow are excluded from this first analysis, because the exponent is intended to describe the transverse valley shape.

The measured beds are represented by

\[
z(x)=z_0+a|x-x_0|^b,
\]

where \(x_0\) and \(z_0\) locate the fitted valley bottom, \(a\) controls the opening of the valley, and \(b\) controls its shape. Values near \(b=1\) describe a comparatively sharp V-shaped valley, whereas \(b=2\) gives a parabola. Across the five retained full cross-sections, the fitted exponent ranges from 1.30 to 1.49, with a mean of 1.36 and a standard deviation of 0.08. The lowest value, and therefore the most V-shaped section in this set, occurs near the glacier tongue.

## Interpreting the free and parabolic fits

Both fits are anchored at the same two points where the bed meets the 2023 surface. Unequal anchor elevations are accommodated by fitting the left and right valley sides separately and joining them at a common valley bottom. With the two anchors fixed, the free fit optimizes the horizontal position and elevation of the valley bottom, \(x_0\) and \(z_0\), together with \(b\); the scaling of each side is determined by its anchor. The parabolic comparison fixes \(b=2\) and optimizes only \(x_0\) and \(z_0\). Both minimize squared vertical residuals over the entire observed cross-section.

RMSE summarizes the average vertical mismatch along the full profile; it does not specifically measure the error at the valley bottom. Two curves can have comparable RMSE values while implying substantially different maximum ice thicknesses.

## Partial GPR cross-sections

We then return to the partial airborne GPR profiles, for which bed reflections are visible on the valley flanks but cannot be identified through the deepest part of the glacier. These profiles do not directly determine the maximum ice thickness. We compare two reconstructions:

1. A transverse power-law curve with \(b=1\). Its remaining parameters are fitted only to the observed GPR bed reflections.
2. A free power-law fit to the Grab et al. (2021) interpolated bed profile.

The Grab bed should not be interpreted as an independent direct observation beneath the GPR gap. It is a model- and data-constrained reconstruction. The GlaTE approach described by [Langhammer et al. (2019)](https://tc.copernicus.org/articles/13/2189/2019/) combines absolute ice-thickness information from GPR with gradients supplied by a glaciological mass-conservation model. In that model component, apparent mass balance and glacier geometry are used to estimate ice discharge through elevation-zone boundaries; ice-flow physics then relates discharge, surface slope, and basal shear stress to thickness. The joint inversion uses measured GPR reflections to constrain absolute thickness where observations exist, while the mass-conservation information regularizes the thickness distribution between and beyond those observations. The resulting Grab profile is thus physically informed, but its deepest section remains dependent on modeling assumptions, regularization, input geometry, mass balance, and ice-flow parameters. In particular, we suspect the smoothing step to underestimate the maximum thickness.

For the partial cross-sections, the \(b=1\) reconstruction and the fit to the Grab bed can diverge substantially in the central, deepest part of the valley. Their RMSE values may nevertheless be similar because the measured GPR constraints occur mainly on the flanks, where both shapes can fit reasonably well. 

## Can gravimetry discriminate between the shapes and better constrain the maximum ice thickness ?

Gravimetry offers a potentially valuable independent constraint because replacing dense bedrock with lower-density ice produces a measurable gravity anomaly. A deeper or broader ice-filled trough contains more low-density material and should generate a different anomaly from a shallower or more V-shaped trough. Forward gravity modeling can therefore test whether candidate geometries are compatible with the observed gravity field.

To remain consistent with a first-order gravimetry approach, we propose inverting only one scalar shape parameter: the exponent **b**. The same trial value is imposed in every partial transverse GPR cross-section. At each transect, the power-law shape is inserted only across the interval where no reliable bed reflection exists. All observed GPR bed elevations on the flanks remain hard constraints (within their uncertainties, typically 5m) and never move. Unequal elevations at the two gap boundaries are accommodated by scaling the left and right sides separately while retaining the same exponent, so the gap fill connects both observed boundaries exactly without requiring a globally symmetric section.

For each trial b, the completed cross-sections are supplied to GlaTE together with all original GPR observations. GlaTE then recomputes the distributed three-dimensional ice-thickness field using the GPR and glaciological modeling constraints. Maximum depth is not incremented or independently inverted at each profile; it emerges from the GlaTE solution associated with the trial exponent. This deliberately simple formulation avoids introducing one depth parameter per transect.

Unfortunately, there is no GPR bed reflection at the two gravimetry profiles. We therefore propose treating each of these cross-sections as part of the parameterized inversion: the section is anchored at the two glacier margins, assigned the same trial exponent \(b\) as the partial GPR sections, and allowed to obtain its maximum depth from the coupled GlaTE–gravity solution. This retains the deliberately low-dimensional inversion and lets gravimetry inform the cross-sections where it has the greatest sensitivity. As a sensitivity test, the two sections can instead be left unconstrained within GlaTE, so their shapes are determined mainly by the glaciological model and spatial smoothing. Comparing the gravity misfit from these two treatments will show how strongly the result depends on prescribing the common cross-sectional shape. Note that we aim to acquire low-frequency GPR measurements at these locations to add direct constraints and validate the approach.

Each updated ice-thickness (or bed-elevation) model is passed to the gravity forward model and evaluated against the observations. We define the gravity misfit as the normalized root-mean-square residual, \(\Phi_g=[N^{-1}\sum_i((g_i^{\mathrm{obs}}-g_i^{\mathrm{pred}})/\sigma_i)^2]^{1/2}\), where \(\sigma_i\) combines measurement uncertainty with the uncertainty assigned to elevation, density contrast, and other forward-model terms. This weighting prevents stations with larger uncertainty from controlling the inversion and permits comparison between trials. The preferred model minimizes \(\Phi_g\); a value near or below one indicates that the residuals are, on average, consistent with their assigned uncertainties. The inversion is stopped when further changes in \(b\) do not materially reduce \(\Phi_g\), or when the residuals fall within the gravimetry uncertainty. The final product is a GPR-, gravimetry-, and glaciological-model-constrained bed. In this exercise, gravity alone will not guarantee a unique bed solution, but it can provide an additional constraint on the GlaTE inversion where direct GPR anchors are absent. Because the inferred bed will depend on the gravity observations and forward-model assumptions, sensitivity tests should quantify the effects of uncertain parameters, such as ice–bedrock density contrast, instrument standard deviation, and elevation. These tests are required to propagate uncertainty into the final bed-elevation product. 

### Active and inactive inversion domains

The adjustment of ice thickness is restricted to an **active domain** along the main Aletsch trunk from Konkordiaplatz to the P21 transect, which is the downstream limit. Furetherdown, the GPR offer bedrock resolved on full cross section and the ice thickness is assumed to be accurate enough (i.e. not possibly better constrained by gravimetry). 

Everything else is the **inactive domain** and the distributed ice thickness values do not move: the bed upstream of Konkordiaplatz, downstream of P21, and in neighboring branches, is taken from the Grab et al (2021) ice thickness product.



### Surface velocities

Observed surface velocities could first be used as an independent validation dataset. For each candidate thickness field, a forward ice-flow calculation would predict surface velocity, and the residuals could reveal locations where the inferred geometry is dynamically inconsistent. In particular, a model that systematically underpredicts velocity in the fast-flowing trunk may indicate underestimated ice thickness, although velocity residuals are not uniquely attributable to thickness because they also depend on basal sliding, ice rheology, surface slope, and mass balance.

Surface velocity could subsequently be incorporated as an additional GlaTE constraint through a separately weighted velocity-misfit term. This should be tested independently from the gravity inversion before combining both datasets. Agreement between velocity- and gravity-supported thickness changes would provide stronger evidence for a deeper bed; disagreement would help expose sensitivity to the glaciological or gravity-model assumptions. Gravimetry would remain useful because it responds directly to the subsurface mass deficit, whereas velocity provides an indirect and potentially non-unique dynamical constraint.

### Validation

We plan more GPR full cross section using low frequency (6 MHz) to validate ice thickness at the gravimetry profile, and in the main glacier trunk (in between the two gravimetry profiles).


## References

- [Langhammer, L., Grab, M., Bauder, A., and Maurer, H. (2019), *Glacier thickness estimations of alpine glaciers using data and modeling constraints*, The Cryosphere, 13, 2189–2202.](https://doi.org/10.5194/tc-13-2189-2019)
- [Grab et al. (2021), *Ice thickness distribution of all Swiss glaciers based on extended ground-penetrating radar data and glaciological modeling*.](https://www.dora.lib4ri.ch/wsl/islandora/object/wsl%3A27301/datastream/PDF/Grab-2021-Ice_thickness_distribution_of_all-%28published_version%29.pdf)
