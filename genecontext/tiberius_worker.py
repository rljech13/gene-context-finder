"""Isolated CPU entry point, executable with the predictor's own Python."""
import argparse
import hashlib
import importlib.metadata as metadata
import json
import os
import sys


def probe():
    packages = {name: metadata.version(name) for name in
                ("tiberius", "tensorflow", "bricks2marble", "hidten", "biopython", "numpy", "keras")}
    files = {str(p): hashlib.sha256(p.locate().read_bytes()).hexdigest()
             for p in metadata.files("tiberius") if str(p).endswith(".py")}
    return dict(python=sys.version, executable=sys.executable, packages=packages,
                tiberius_source_sha256=files)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--probe", action="store_true")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--model")
    parser.add_argument("--fasta")
    parser.add_argument("--out")
    parser.add_argument("--batch-size", type=int, default=1)
    args = parser.parse_args()
    if args.probe:
        print(json.dumps(probe()))
        return
    if not args.model or not args.fasta or not args.out or min(args.threads, args.batch_size) < 1:
        parser.error("model, fasta, out and positive thread/batch counts are required")
    os.environ.update(CUDA_VISIBLE_DEVICES="-1", TF_NUM_INTRAOP_THREADS=str(args.threads),
                      TF_NUM_INTEROP_THREADS="1", OMP_NUM_THREADS=str(args.threads),
                      OPENBLAS_NUM_THREADS=str(args.threads), TF_CPP_MIN_LOG_LEVEL="1")
    import tensorflow as tf
    tf.config.set_visible_devices([], "GPU")
    tf.config.threading.set_intra_op_parallelism_threads(args.threads)
    tf.config.threading.set_inter_op_parallelism_threads(1)
    from tiberius.main import run_tiberius, clamp_seq_len_to_genome
    from tiberius.tiberius_args import parseCmd
    sys.argv = [__file__, "--genome", args.fasta, "--model", args.model,
                "--no_softmasking", "--seq_len", "400050", "--batch_size", str(args.batch_size),
                "--out", args.out, "--codingseq", args.out + ".cds.fna",
                "--protseq", args.out + ".proteins.faa"]
    settings = parseCmd()
    effective = clamp_seq_len_to_genome(args.fasta, settings.seq_len, settings.min_genome_seqlen)
    with open(args.fasta) as handle:
        input_windows = sum(line.startswith(">") for line in handle)
    print("GENE_CONTEXT_CPU=" + json.dumps(dict(tensorflow=tf.__version__,
          visible_devices=[str(device) for device in tf.config.get_visible_devices()],
          intraop_threads=args.threads, interop_threads=1, input_windows=input_windows,
          requested_seq_len=settings.seq_len, effective_seq_len=effective)), flush=True)
    run_tiberius(settings)


if __name__ == "__main__":
    main()
