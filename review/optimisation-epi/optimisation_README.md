# Documentation for optimisation files

## Description of files

`single_run.R`

Description:

Runs SEIR for 360 days - with primary emphasis on the 2 vaccines v2 model. More details on model specification and parameters can be found in `review.pdf`. Calls numerical integration code from `setup/model_2vaccine.R`.

Functions:
`sr` : outputs numerical integration results from `setup/model_2vaccine.R` given initial conditions and model parameters.

Depends on:
- `setup/model_2vaccine.R`

`objective-functions.R`

Description:
Defines objective functions for non-linear optimization problem given the `single_run.R` defined above - example objectives include total deaths at day 360. 

Functions:
`rescale_rcs` - turns proportional outputs into population outputs
`phi_x` - dose response function
`model_fd_dynamic` - dynamic model (i.e. different vaccine allocations to different age groups)
`model_fd_static` - static model (i.e. uniform vaccine allocations to different age groups)

`nlopt-analysis.R`

Description:
Runs optimization scripts (optimizing over allocation by age group) over the objective functions `objective-functions.R`. 

Outputs:
- allocations to `results/wip-nl-solutions7*.Rdata`

Depends on:
- `project_setup.R`
- `optimisation-epi/objective-functions.R`

`visualise-results.R`


