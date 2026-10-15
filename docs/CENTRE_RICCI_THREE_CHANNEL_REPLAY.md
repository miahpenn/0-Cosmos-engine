# Centre-Ricci Lambda three-channel replay

**Purpose:** separate three derivative channels in the initialized/evolved-state centre-Ricci investigation without starting another production-kernel run.

The offline replayer is \`tools/decompose_saved_centre_ricci_profiles.py\`. It consumes a saved stage-budget artifact ZIP, verifies the outer ZIP digest and the SHA-256 manifest inside it, checks the artifact's source/vendor pins against the current audit, then writes a JSON summary and an NPZ of the decomposed profiles. It never initializes or advances the production kernel.

## Reproducing the recorded N=80 stage capture

Use the checksum-verified artifact from [stage-budget run 38692095225](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38692095225), artifact ID \`11952561767\`. Its outer ZIP SHA-256 is:

\`7f64d7fd15347c25924f83381133a6df4347c692239299ac1f58c2827f43ca0c\`

Run:

\`\`\`bash
python tools/decompose_saved_centre_ricci_profiles.py \\
  --artifact /path/to/ricci-stage-n80.zip \\
  --expected-zip-sha256 7f64d7fd15347c25924f83381133a6df4347c692239299ac1f58c2827f43ca0c \\
  --output results/current-main-centre-ricci-three-channel/n80.json
\`\`\`

The expected production-source anchor is \`285a2ee0c1a3b97cddeb4519d8e749109b115aab\`; the vendor pin is \`d6052d605673ce9d82cdc99b0df79855fa2ea215\`. The source archive manifest, source pin and vendor pin are mandatory checks, not annotations.

## Channel definition

For the same sampled metric, let \`Lambda_evolved\` be the saved independently evolved connection and \`Lambda_metric\` be reconstructed from the sampled metric with the grid's \`D1\` operator. Let \`D4\` mean the pinned grid's fourth-order derivative operator. The report keeps these terms separate:

- \`G0\`: the vendor TF-Ricci gap using \`D4(Lambda_evolved)\`.
- \`C = (2/3) X^2 [D4(Lambda_metric) - D4(Lambda_evolved)]\`: the sampled connection-residual / D4 channel.
- \`S = (2/3) X^2 [partial_r(Lambda_metric)|D1/D2 - D4(Lambda_metric)]\`: the product/chain-rule stencil channel.
- \`Gchain = G0 + C + S\`: must reproduce the previously saved chain-rule counterfactual gap to round-off.

The regression first exercises this exact recomposition on smooth, regular-centre manufactured profiles, including a deliberately off-constraint odd \`Lambda\` contribution. A second regression round-trips a synthetic artifact through the offline reader and asserts the replay executes zero evolution steps. The N=80 artifact can then be reprocessed without running the one-step capture again.

## Interpretation boundary

A small \`Gchain\` or tiny identity-closure residual confirms only the internal arithmetic and channel attribution on these samples. It does not prove that the BSSN Ricci representation is physically wrong, authorize changing the evolved production \`Lambda\`, prove long-horizon stability, or constitute scientific admission. Keep production/vendor equations unchanged pending independent scientific review and the missing first-principles source/regularity and convergence gates.
