### RVF.Vaxx

Rift Valley Fever is modeled as a two-host/two-vector system, with disease
transmission occurring between humans, a mammalian livestock species, Culex
mosquitos, and Aedes mosquitos.

Humans and livestock progress through four disease stages of: susceptible,
exposed, infected, and recovered. Waning immunity from natural infection is
assumed. Vaccination is possible and occurs without
regard to disease status. There is a two-week delay before vaccination is
protective and waning immunity from vaccination is also assumed, requiring
booster vaccinations every year for livestock and every five years for humans.

Mosquitos of both species include adult disease states of susceptible and
infectious. Both disease states lay eggs with a similar rate. Infectious vectors
convey infection to eggs probabilistically.

&nbsp;

#### Package Contents:

`RVF.Vaxx` currently contains two functions:

* `rvf_vaxx()` models RVF transmission as a series of 20
deterministic ordinary differential equations representing movement between
disease state compartments.

* `experiment_run()` allows for tracked simple experiments to manipulate initial
states and rate of change parameters, with an auto-generated diagnostic plot to
visualize effects of experimentation.

&nbsp;

#### Model Assumptions:

1.  All human/mammalian births are susceptible; no conferred immunity.

2.  No testing before vaccination; all exposure categories vaccinated at same rate.

3.  Vaccinating exposed humans/mammals prevents conversion to infection.

4.  Waning immunity from recovery and vaccination.

5.  Adjusted birth rate for mammalian host assumes differential births due to maternal infection.

6.  Vector hatch rate is asssumed to be independent of disease status.
