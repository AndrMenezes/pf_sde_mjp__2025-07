
Use of SMC-like algorithms to try to improve inference and filtering on high
(or, even moderate) dimension Markov jump processes.

## Some initial references:
There are some really messy papers on using SMC for filtering on SPDEs,
but the process is so high-dimensional that they don't even try to then do
inference
(I don't recommend reading them as a first pass - there are also many typos in
one of them that just add to the confusion).
- https://arxiv.org/abs/1710.04586
- https://arxiv.org/abs/1907.11884

There's also another paper that sets out a general framework; but it is very abstract and too generic.
https://warwick.ac.uk/fac/sci/statistics/staff/academic-research/johansen/publications/btpf.pdf


> Chris is putting in a grant proposal on new methods for MJPs, but even if it comes off that wouldn't start until towards the end of 2026.

To look at adapting the overall idea in the above to permit filtering and inference for moderate-dimensional MJPs.
I already have some specific thoughts on the changes that would be needed, so the project could hit the ground running.

## Project tasks

- Coding up a simple version of a particle MCMC scheme for inference on a
moderate-dimensional (time-discretised) SDE.

- Developing tools to do this for an MJP (simulated using the tau-leap scheme),

- Seeing how high we can push the dimension e.g. for a set of partially observed,
coupled SEIR models, where you have fairly accurate observations of the number of
new infecteds in each age category or health/spatial region.

- There may be some interesting work to do on tuning the correlation parameter
 for the proposed move (I'm guessing that choosing the tempering according to
 the standard adaptive scheme will probably work fine, so we probably won't
 have to worry about that).

## Comments

Implicit, I need to understand:

- Particle filter (bootstrap one).

- Sthocastic diferential equation.

- Particle MCMC
  - Here I need to ask if we do MCMC only for parameter inference of for the filtering as well:
    - Awnser: only for parameter inference.

- Tempering, as for the MJPs and SDE they used tempering. And Chris is also asking for do that.


Given particle MCMC is completely new to you, we should take a step back from my initial suggestion.  It would be sensible to start off with coding up a bootstrap particle filter for a simple system (e.g. try the linear Gaussian model, where you can work out the true likelihood via the Kalman filter and, hence, check your particle filter code is correct) and then using this to perform particle MCMC for inference on the parameters of the system.

- Use this library: https://www.boost.org/
