## R CMD check results

0 errors | 0 warnings | 0 notes

(If submitted before the server-side `dsRctFed` package is on CRAN, the
check additionally produces one NOTE: "Suggests or Enhances not in
mainstream repositories: dsRctFed". See below.)

## Companion server-side package

`dsRctFedClient` is the client-side half of a two-package 'DataSHIELD'
pair. The companion server-side package `dsRctFed` is submitted first;
`dsRctFedClient` is submitted immediately after `dsRctFed` is accepted.

`dsRctFed` is listed under `Suggests:` because the shared `getting_started`
vignette and some tests exercise the end-to-end federated workflow against
local `DSLite` servers that load `dsRctFed`. Those chunks and tests are
guarded with `requireNamespace()` and are skipped automatically when
`dsRctFed` is not installed.

## Test environments

* local Ubuntu, R 4.6.0
* win-builder (devel and release)

## Reverse dependencies

None (new package).
