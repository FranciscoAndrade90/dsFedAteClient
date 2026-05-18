# dsRctFedClient 0.0.0.9000

* Renamed `connections` argument to `datasources` in `ds.EstimateAte()`.
* Added `datashield.connections_find()` fallback when `datasources = NULL`.
* Added `.check_datasources()` to validate the datasources argument.
* Added `.check_args()` and `.check_string()` to validate character arguments.
* Renamed internal helper `ds.localAte()` to `.localAte()` and marked it internal.
* Removed `suppressWarnings()` wrappers from `ds.glm()` calls.
* Added `cli` and `methods` under `Imports` in `DESCRIPTION`.
* Added unit tests for `ds.EstimateAte()` and all input validation helpers.
