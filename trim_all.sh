#!/usr/bin/env bash
# Trimmomatic on all runs in raw_reads/, 4 runs at a time (2 threads each).
# Skips runs that already finished, so it is safe to rerun after a restart.
# Per-run log: trimmed/<run>/<run>.trim.log   Summary: trimmed/trim_summary.tsv
# Usage (inside typhi_amr env):  nohup bash trim_all.sh > trim.log 2>&1 &

set -u
cd "$(dirname "$0")"
AD=$CONDA_PREFIX/share/trimmomatic-0.40-0/adapters/TruSeq3-PE-2.fa
[ -f "$AD" ] || { echo "adapter file not found: $AD (is typhi_amr active?)"; exit 1; }
mkdir -p trimmed

trim_one() {
    r=$1; o=trimmed/$r
    if [ -f $o/$r.done ]; then echo "skip $r"; return; fi
    mkdir -p $o
    if trimmomatic PE -threads 2 -phred33 \
        raw_reads/$r/${r}_1.fastq.gz raw_reads/$r/${r}_2.fastq.gz \
        $o/${r}_1P.fastq.gz $o/${r}_1U.fastq.gz $o/${r}_2P.fastq.gz $o/${r}_2U.fastq.gz \
        ILLUMINACLIP:$AD:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50 \
        > $o/$r.trim.log 2>&1; then
        touch $o/$r.done; echo "ok $r"
    else
        echo "FAILED $r (see $o/$r.trim.log)"
    fi
}
export -f trim_one; export AD

ls raw_reads | xargs -P 4 -I{} bash -c 'trim_one {}'

# summary table from the "Input Read Pairs" line of each log
echo -e "run\tinput_pairs\tboth_surviving\tpct_both\tdropped_pct" > trimmed/trim_summary.tsv
for f in trimmed/*/*.trim.log; do
    r=$(basename $f .trim.log)
    grep '^Input Read Pairs' $f | awk -v r=$r '{gsub(/[()%]/,""); print r"\t"$4"\t"$7"\t"$8"\t"$21}'
done >> trimmed/trim_summary.tsv
echo "done: $(ls trimmed/*/*.done | wc -l)/$(ls raw_reads | wc -l) runs trimmed"
