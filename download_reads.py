#!/usr/bin/env python3
"""Download all paired FASTQs for PRJEB6352 from the ENA manifest, in parallel."""
import csv
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

MANIFEST = "prjeb6352_final_manifest.csv"
OUTDIR = "raw_reads"
PARALLEL = 8

def build_jobs():
    jobs = []
    with open(MANIFEST, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            run = row["run_accession"]
            ftp_field = row["fastq_ftp"]
            if not ftp_field:
                print(f"WARNING: no fastq_ftp for {run}", file=sys.stderr)
                continue
            run_dir = os.path.join(OUTDIR, run)
            os.makedirs(run_dir, exist_ok=True)
            for url in ftp_field.split(";"):
                url = url.strip()
                if not url:
                    continue
                if url.startswith("ftp://"):
                    url = "http://" + url[6:]
                elif not url.startswith("http"):
                    url = "http://" + url
                dest = os.path.join(run_dir, os.path.basename(url))
                jobs.append((url, dest))
    return jobs

def fetch(job):
    url, dest = job
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return (dest, "skipped-exists")
    tmp = dest + ".part"
    wget_args = [
        "wget",
        "--tries=5",
        "--waitretry=10",
        "--timeout=60",
        "--read-timeout=60",
        "--retry-on-http-error=403,429,500,502,503,504",
        "--user-agent=typhi-study-downloader/1.0",
    ]
    if os.path.exists(tmp):
        wget_args.append("--continue")
    wget_args.extend(["-O", tmp, url])
    result = subprocess.run(
        wget_args,
        capture_output=True, text=True,
    )
    if result.returncode == 0:
        os.rename(tmp, dest)
        return (dest, "ok")
    else:
        return (dest, f"FAILED: {result.stderr.strip()[:200]}")

if __name__ == "__main__":
    jobs = build_jobs()
    print(f"Total files to download: {len(jobs)}")
    ok, skipped, failed = 0, 0, []
    with ThreadPoolExecutor(max_workers=PARALLEL) as pool:
        futures = {pool.submit(fetch, job): job for job in jobs}
        done = 0
        for fut in as_completed(futures):
            dest, status = fut.result()
            done += 1
            if status == "ok":
                ok += 1
            elif status == "skipped-exists":
                skipped += 1
            else:
                failed.append((dest, status))
            print(
                f"[{done}/{len(jobs)}] {status}: {dest} "
                f"(ok={ok} skipped={skipped} failed={len(failed)})",
                flush=True,
            )
    print("\n=== Summary ===")
    print(f"OK: {ok}  Skipped (already present): {skipped}  Failed: {len(failed)}")
    if failed:
        print("\nFailed downloads:")
        for dest, status in failed:
            print(f"  {dest}: {status}")
        sys.exit(1)
