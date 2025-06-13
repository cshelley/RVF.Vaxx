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

![Aly_Figure1](https://github.com/user-attachments/assets/8aabfde6-842a-48bc-b63c-06c14f3e03ac)


&nbsp;

#### Package Contents:

`RVF.Vaxx` currently contains three functions:

* `rvf_vaxx()` models RVF transmission as a series of 20
deterministic ordinary differential equations representing movement between
disease state compartments.

* `experiment_run()` allows for tracked simple experiments to manipulate initial
states and rate of change parameters, with an auto-generated diagnostic plot to
visualize effects of experimentation.

* `rvf_sensitivity()` is a simplified version of the full model prepared for
sensitivity analysis by removing non-actionable compartments and simplifying
parameters.

&nbsp;

#### Full Model Assumptions:

1.  All human/mammalian births are susceptible; no conferred immunity.

2.  No testing before vaccination; all exposure categories vaccinated at same rate.

3.  Vaccinated individuals move to vaccination compartment, stopping further
disease progression.  

4.  Waning immunity from recovery and vaccination.

5.  Adjusted birth rate for mammalian host assumes differential births due to maternal infection.

6.  Vector hatch rate is asssumed to be independent of disease status.

&nbsp;

![Aly_Figure2](https://github.com/user-attachments/assets/cb593e54-52f7-43b7-aa3b-9060f4aae92e)

#### Sensitivity Model Assumptions: 

1. No births/deaths in human/mammalian hosts.

2. Mosquito egg laying/death rates are not affected by disease status.
