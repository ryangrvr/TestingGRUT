# Implementation attempts

The first short-check import failed because the new module name `common` collided with a legacy reheating dependency. It was renamed `pt_common` before any new calculation completed. The first spectral analysis likewise found a legacy `production` module; the new streaming driver now has the unique name `pt_production`. Those import failures produced no scientific results and changed no model inputs.

A missing closing parenthesis in the response driver was caught before execution and corrected. The working response uses the same scalar Hamiltonian and its exact discrete tangent. No scientific threshold was adjusted in response to these implementation failures.

The two original-band production workers were already running when the driver filename was changed. Their loaded stepping and observation code was unchanged; the remaining infrared runs use the renamed driver. Both names refer to the same numerical implementation, and the reproduction commands use the final name.
