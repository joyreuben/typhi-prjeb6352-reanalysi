#!/usr/bin/env bash
# Typhi Mykrobe (GenoTyphi genotype + QRDR mutations + AMR genes + plasmids) on the 159 Typhi runs,
# from trimmed reads. 2 runs at a time x 2 threads. Skips runs that already have a JSON, so safe to rerun.
# Output: mykrobe/<run>.json   Summary table: mykrobe/mykrobe_out_predictResults.tsv
# Usage:  nohup bash mykrobe_all.sh > mykrobe.log 2>&1 &

set -u
cd "$(dirname "$0")"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mykrobe
mkdir -p mykrobe

mykrobe_one() {
    r=$1
    if [ -s mykrobe/$r.json ]; then echo "skip $r"; return; fi
    if mykrobe predict --sample $r --species typhi --format json --threads 2 \
        --tmp mykrobe/tmp_$r --output mykrobe/$r.json.part \
        --seq trimmed/$r/${r}_1P.fastq.gz trimmed/$r/${r}_2P.fastq.gz > mykrobe/$r.stdout 2>&1 \
        && [ -s mykrobe/$r.json.part ]; then
        mv mykrobe/$r.json.part mykrobe/$r.json; rm -rf mykrobe/tmp_$r mykrobe/$r.stdout; echo "ok $r"
    else
        echo "FAILED $r (see mykrobe/$r.stdout)"
    fi
}
export -f mykrobe_one

xargs -P 2 -I{} bash -c 'mykrobe_one {}' < results/typhi_runs.txt
echo "done: $(ls mykrobe/*.json | wc -l)/$(wc -l < results/typhi_runs.txt) runs"
conda activate typhi_amr
python genotyphi/typhimykrobe/parse_typhi_mykrobe.py --jsons mykrobe/*.json --prefix mykrobe/mykrobe_out
