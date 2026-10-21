# Claude's audit of Codex's 3D charged-exit cross-check (a69bd86)

Start with [REVIEW.md](REVIEW.md).
- [BRIEF.md](BRIEF.md) is the shared brief given to the audit agents.
- `audit/<track>/` holds each track's NOTES.md, code and results; `audit/v_<track>/` holds the adversarial verifier
  of that track.

Tracks:
- **rerun** reproduces Codex's run table and reruns its frozen C++ solver, including the two skipped large-box seeds.
- **preparation** explains the one-vacuum sparse seed.
- **estimators** checks wall-speed and Lorentz weighting.
- **wbound** works out what can be proven about the hidden equation of state after parity.
- **reach** sets out the empirical channels of a gravity-only hidden sector.

The corrections this audit accepted are applied in `../grut_issue13_wall_history_20261020/REPORT.md` §11.
