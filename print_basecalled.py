import argparse
import os
import subprocess
import sys
import tempfile
import time

import pod5


def parse_args():
    parser = argparse.ArgumentParser(
        description="Watch a growing .pod5 file and basecall each new read as it appears."
    )
    parser.add_argument("path_to_dorado_models", help="Path to the Dorado basecalling model directory")
    parser.add_argument("pod5_file", help="Path to the .pod5 file being written by MinKNOW")
    parser.add_argument(
        "--poll-interval",
        type=float,
        default=2.0,
        help="Seconds to wait between checks for new reads (default: 2.0)",
    )
    return parser.parse_args()


def run(command, error_message):
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"{error_message}\nstderr: {result.stderr}")


def basecall_read(read, path_to_dorado_models, work_dir):
    tmp_pod5 = os.path.join(work_dir, "tmp.pod5")
    tmp_bam = os.path.join(work_dir, "tmp.bam")
    tmp_fasta = os.path.join(work_dir, "tmp.fasta")

    # write the single new read to its own pod5 file
    with pod5.Writer(tmp_pod5) as writer:
        writer.add_read(read)

    # basecall it
    with open(tmp_bam, "wb") as bam_out:
        result = subprocess.run(
            ["dorado", "basecaller", path_to_dorado_models, tmp_pod5],
            stdout=bam_out,
            stderr=subprocess.PIPE,
            text=False,
        )
    if result.returncode != 0:
        raise RuntimeError(f"Failed to basecall {read.read_id}!\nstderr: {result.stderr.decode(errors='replace')}")

    # convert to fasta
    with open(tmp_fasta, "w") as fasta_out:
        result = subprocess.run(
            ["samtools", "fasta", tmp_bam],
            stdout=fasta_out,
            stderr=subprocess.PIPE,
            text=True,
        )
    if result.returncode != 0:
        raise RuntimeError(f"Failed to convert to fasta {read.read_id}!\nstderr: {result.stderr}")

    with open(tmp_fasta) as fasta_file:
        content = fasta_file.read()
        if content:
            print(content, flush=True)

    # clean up this read's temp files before the next one
    for path in (tmp_pod5, tmp_bam, tmp_fasta):
        try:
            os.remove(path)
        except OSError:
            pass


def main():
    args = parse_args()
    seen_read_ids = set()

    print(f"Watching {args.pod5_file} for new reads... (Ctrl+C to stop)", flush=True)

    with tempfile.TemporaryDirectory(prefix="basecall_") as work_dir:
        try:
            while True:
                with pod5.Reader(args.pod5_file) as reader:
                    # .to_read() must be called while the reader is still open —
                    # ReadRecord objects are lazy references into the Arrow tables,
                    # which are released when the context exits.
                    reads = [(r.read_id, r.to_read()) for r in reader.reads()]

                new_reads = [(rid, r) for rid, r in reads if rid not in seen_read_ids]
                
                for read_id, read in new_reads:
                    basecall_read(read, args.path_to_dorado_models, work_dir)
                    seen_read_ids.add(read_id)

                if not new_reads:
                    total = len(reads)
                    print(f"Čekám na nové čtení... (zatím přečteno {total})", flush=True)
                time.sleep(args.poll_interval)
        except KeyboardInterrupt:
            print("\nStopped.")
            sys.exit(0)


if __name__ == "__main__":
    main()