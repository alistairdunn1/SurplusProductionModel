# SurplusProductionModel Development Roadmap

## Milestone v0.1.1 (Stabilization and Consistency)
Target window: 1-2 weeks

### 1. Unify reference point math across modules [IN PROGRESS]
Effort: M (1-2 days)

Tasks:
- Replace duplicate MSY/BMSY/FMSY math in fitting internals with shared helper.
- Confirm profile and Bayesian derived values use shared helper outputs.
- Document one canonical formula source in code comments.

Acceptance criteria:
- One canonical reference point helper is used by all production paths.
- Existing tests pass with no new failures.

### 2. End-to-end depletion target consistency tests
Effort: M (1 day)

Tasks:
- Add cross-function test: calculate_reference_points vs profile_likelihood vs bayesian_fit for same depletion targets.
- Add status/print coverage for user-defined targets in summary output.

Acceptance criteria:
- Depletion target values are numerically consistent across pathways.

### 3. Profile likelihood warning hardening
Effort: S-M (0.5-1 day)

Tasks:
- Handle invalid profile grid points gracefully (mark non-converged/NA instead of warning spam).
- Reduce NA/NaN evaluation warnings surfaced during test runs.

Acceptance criteria:
- No regression in CI results.
- Warning volume reduced meaningfully in profile tests.

### 4. Documentation and examples sync
Effort: S (0.5 day)

Tasks:
- Align README and Rd docs with current arguments and behavior.
- Add one compact example for depletion targets in profile and Bayesian workflows.

Acceptance criteria:
- Help pages and README reflect current API exactly.

## Milestone v0.2.0 (Reliability and Release Readiness)
Target window: 2-4 weeks

### 5. Bayesian CI coverage on dependency-enabled runner [IN PROGRESS]
Effort: M (1 day)

Tasks:
- Add a CI job (or scheduled workflow) with tmbstan/rstan available.
- Enable Bayesian path tests there.

Acceptance criteria:
- Bayesian tests execute regularly, not only skipped locally.

### 6. User-facing defaults for target reporting
Effort: M (1 day)

Tasks:
- Add optional package-level defaults for biomass_target and baseline.
- Keep function-level override behavior explicit.

Acceptance criteria:
- Users can set defaults once and get consistent output across summary functions.

### 7. Release hygiene and changelog
Effort: S (0.5 day)

Tasks:
- Run full checks and resolve blocking warnings/errors.
- Add NEWS entry for depletion-target feature set and consistency refactor.

Acceptance criteria:
- Package ready for tagged release.

## Immediate Next Action (this session)
- Start item 1 by replacing the fitting-internal reference point calculation with the shared helper.
