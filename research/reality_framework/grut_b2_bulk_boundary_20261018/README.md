# B2-BB actual bulk and boundary kernel

Start with [REPORT.md](REPORT.md). The actual spatial reduction is
author-derived; the original operating gate is unresolved.

- [BOUNDS.md](BOUNDS.md): verifies derivative/monotonicity and source tail;
  preserves the stronger existing high-tail bound.
- [BULK_BOUNDARY.md](BULK_BOUNDARY.md): derives the continued fraction,
  complete boundary-word bound and trace-norm margin.
- [SOURCE_KERNEL.md](SOURCE_KERNEL.md): actual spatial inverse from
  Fourier-log factors, exact source/fuel coupling and remaining theorem.
- [REVIEW_PROBLEM.md](REVIEW_PROBLEM.md): focused real-axis scattering task.
- [EXECUTION_NOTES.md](EXECUTION_NOTES.md): what ran and what did not.

Charter ecff27de1be31a3a383b89ea63dfdee271f91d67; proof/code freeze
ca9ce9a00f834e686a7c68e74883b876517670f7. No historical law or gate changed.
Run `python controls.py`, then `python verify.py`. Exact control agreement
does not settle spectral type or source operation. No follow-on automated
campaign is opened by this packet.
