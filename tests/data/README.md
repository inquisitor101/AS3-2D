# Test data

`grid_AS3/zone_0.as3` through `zone_8.as3` and `param_t1.cfg` are the unchanged
fixtures transferred from the original Skylabs tests. Their checksums are listed
in `../reference_manifest.json`.

The numerical unit tests construct their inputs directly and do not depend on
the working directory or a simulation configuration. This parameter file records
the original fixture setup; it uses earlier configuration names and is not
currently registered as a simulation regression input.
