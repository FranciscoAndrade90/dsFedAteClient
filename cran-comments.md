## R CMD check results

0 errors | 0 warnings | 1 note

The NOTE combines two expected items:

* New submission.
* Suggests or Enhances not in mainstream repositories: dsFedAte.

See below.

## Companion server-side package

`dsFedAteClient` is the analyst-side half of a two-package 'DataSHIELD'
pair. The companion server-side package `dsFedAte` is submitted first;
`dsFedAteClient` is submitted immediately after `dsFedAte` is accepted,
at which point the "Suggests not in mainstream repositories" NOTE
clears.

`dsFedAte` is listed under `Suggests:` because the `federated_rct`
vignette and the test suite exercise the end-to-end federated workflow
against local `DSLite` servers that load `dsFedAte`. Those vignette
chunks and tests are guarded with `requireNamespace()` and are skipped
automatically when `dsFedAte` is not installed, so `R CMD check` on
CRAN's build systems runs cleanly even before `dsFedAte` is published.

## Test environments

* local Ubuntu, R 4.6.0
* win-builder (devel and release)

## Reverse dependencies

None (new package).
