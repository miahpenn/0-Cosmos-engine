#!/bin/bash
# Stage 2 (Jeremiah's go 08:42:29 UTC). Phase arg: inputs | arms | score
set -u
cd /home/user/0-cosmos-engine
VPY=/tmp/claude-0/-home-user-0-cosmos-engine/5354bfc8-ba5f-5b5a-91e5-93af2b317ecf/scratchpad/venv31213/bin/python
OUT=/tmp/claude-0/-home-user-0-cosmos-engine/5354bfc8-ba5f-5b5a-91e5-93af2b317ecf/scratchpad/stage2
export PYTHONPATH=. OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p $OUT
case $1 in
inputs)
  $VPY -m research.control_tests.test2 inputs $OUT/inputs > $OUT/inputs.log 2>&1; echo "inputs exit $?"
  python3 -m research.control_tests.mp_readback $OUT > $OUT/readback_inputs.log 2>&1; echo "readback exit $?";;
arms)
  for a in T2-3 T2-4 T2-6 T2-E; do
    ( $VPY -m research.control_tests.test2 arm $a $OUT > $OUT/$a.log 2>&1; echo "$a $?" >> $OUT/exit_codes.txt ) &
  done; wait; cat $OUT/exit_codes.txt;;
score)
  $VPY -m research.control_tests.test2 cdot_floor results/research/control_tests_stage1/arms $OUT/cdot_floor.json > $OUT/cdot_floor.log 2>&1; echo "cdot_floor exit $?"
  python3 -m research.control_tests.mp_readback $OUT > $OUT/readback.log 2>&1; echo "readback exit $?"
  python3 -m research.control_tests.score2 results/research/control_tests_stage1 $OUT $OUT/cdot_floor.json > $OUT/score2.log 2>&1; echo "score2 exit $?";;
esac
