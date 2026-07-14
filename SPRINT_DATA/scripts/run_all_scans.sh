#!/bin/zsh
cd /Users/saned/dev/SPRINT
export OMP_NUM_THREADS=$(sysctl -n hw.ncpu)
for d in BioGRID HPRD IntAct MINT innate_Exp Innate_Man; do
  echo ">>> START vs-all $d $(date)"
  ./bin/predict_interactions -p human/human_plus_nsp13.fasta -h human/hsp_plus_nsp13 \
    -tr human/PPI_dataset/$d/train.pos.16.txt -pos human/test_nsp13_vs_all.txt \
    -o human/scan_${d}.txt 2>&1 | tail -1
  echo ">>> DONE vs-all $d exit=$? $(date)"
done
echo ">>> ALL SCANS DONE"
