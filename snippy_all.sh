#!/usr/bin/env bash
# Snippy vs CT18 chromosome on all trimmed runs, 2 runs at a time (4 CPUs, 3 GB RAM each).
# Skips runs that already finished, so it is safe to rerun after a restart.
# Output: snippy/<run>/   Per-run log: snippy/<run>/<run>.log
# Usage (inside typhi_amr env):  nohup bash snippy_all.sh > snippy.log 2>&1 &

set -u
cd "$(dirname "$0")"
[ -f ref/CT18.gbk ] || { echo "ref/CT18.gbk missing"; exit 1; }
command -v snippy >/dev/null || { echo "snippy not found (is typhi_amr active?)"; exit 1; }
mkdir -p snippy

snippy_one() {
    r=$1; o=snippy/$r
    if [ -f $o/$r.done ]; then echo "skip $r"; return; fi
    if snippy --force --cpus 4 --ram 3 --ref ref/CT18.gbk \
        --R1 trimmed/$r/${r}_1P.fastq.gz --R2 trimmed/$r/${r}_2P.fastq.gz \
        --outdir $o --prefix $r > snippy/$r.stdout 2>&1 && [ -s $o/$r.tab ]; then
        touch $o/$r.done; rm -f snippy/$r.stdout; echo "ok $r"
    else
        echo "FAILED $r (see snippy/$r.stdout)"
    fi
}
export -f snippy_one

ls trimmed | grep -v '\.tsv$' | xargs -P 2 -I{} bash -c 'snippy_one {}'
echo "done: $(ls snippy/*/*.done | wc -l)/$(ls trimmed | grep -vc '\.tsv$') runs"
