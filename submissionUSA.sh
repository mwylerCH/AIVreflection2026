#!/bin/bash
#SBATCH --job-name="USgeno"
#SBATCH --time=96:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --mail-user=michele.wyler@unibe.ch

# module loading
TEMPDIR="$HOME/temp_AIVms"
export PATH="$PATH:$HOME"

ls $TEMPDIR/singleGeno/*.fa | while read IDX; do
  NAME=`basename $IDX .fa`
  genoflu.py -f $IDX -n $TEMPDIR/genoflu_out/USclass_$NAME &&
  mv $IDX ${IDX}_done
done

