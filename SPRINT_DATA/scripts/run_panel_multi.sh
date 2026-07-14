#!/bin/zsh
cd /Users/saned/dev/SPRINT
export OMP_NUM_THREADS=$(sysctl -n hw.ncpu)
for d in HPRD IntAct; do
  echo ">>> START $d $(date)"
  ./bin/predict_interactions \
    -p human/human_plus_nsp13.fasta \
    -h human/hsp_plus_nsp13 \
    -tr human/PPI_dataset/$d/train.pos.16.txt \
    -pos human/test_nsp13_panel.txt \
    -o human/result_panel_$d.txt 2>&1 | tail -2
  echo ">>> DONE $d exit=$? $(date)"
done
echo ">>> ALL DONE"
