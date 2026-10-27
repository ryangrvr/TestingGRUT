# Environment for track pt_relax (Claude). Codex code runs only from the private copy $RF (git archive of 5f5513b).
export SP=/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad
export RF=$SP/review/pt_relax/rf
export GRUT_KINETIC_INPUT_ROOT=$SP/pinned/research/reality_framework/grut_formation_review_20261025
export GRUT_KINETIC_WORKDIR=$SP/kwork GRUT_ONSET_WORKDIR=$SP/owork
export GRUT_TIMING_WORKDIR=$SP/review/pt_relax/twork GRUT_PARTICLE_WORKDIR=$SP/review/pt_relax/pwork
export GRUT_PARTICLE_WORKERS=1 GRUT_PARTICLE_THREADS=2 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2
export PK=$RF/grut_particle_transport_20261027
export RUNDIR=$SP/review/pt_relax/run
