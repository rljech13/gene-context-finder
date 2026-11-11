"""Isolated CPU entry point, executable with the predictor's own Python."""
import argparse
from contextlib import redirect_stdout, redirect_stderr
import time
import traceback
import hashlib
import importlib.metadata as metadata
import json
import logging
import os
import sys


def probe():
    packages = {name: metadata.version(name) for name in
                ("tiberius", "tensorflow", "bricks2marble", "hidten", "biopython", "numpy", "keras")}
    files = {str(p): hashlib.sha256(p.locate().read_bytes()).hexdigest()
             for p in metadata.files("tiberius") if str(p).endswith(".py")}
    return dict(python=sys.version, executable=sys.executable, packages=packages,
                tiberius_source_sha256=files)


def predict(args, state=None):
    try:
        return _predict(args, state)
    finally:
        # bricks2marble replaces root handlers only when annotation starts.
        # The next window's model-loading messages must not append to an
        # already sealed previous window's log.
        logger = logging.getLogger()
        for handler in list(logger.handlers):
            if isinstance(handler, logging.FileHandler):
                handler.flush(); handler.close(); logger.removeHandler(handler)


def _predict(args, state=None):
    import tensorflow as tf
    from tiberius.main import run_tiberius, clamp_seq_len_to_genome
    from tiberius.tiberius_args import parseCmd
    if state and state.get('path') != args.model:
        state.clear()
        tf.keras.backend.clear_session()
    previous_argv = sys.argv
    sys.argv = [__file__, "--genome", args.fasta, "--model", args.model,
                "--no_softmasking", "--seq_len", "400050", "--batch_size", str(args.batch_size),
                "--out", args.out, "--codingseq", args.out + ".cds.fna",
                "--protseq", args.out + ".proteins.faa"]
    try:
        settings = parseCmd()
    finally:
        sys.argv = previous_argv
    effective = clamp_seq_len_to_genome(args.fasta, settings.seq_len, settings.min_genome_seqlen)
    with open(args.fasta) as handle:
        input_windows = sum(line.startswith(">") for line in handle)
    print("GENE_CONTEXT_CPU=" + json.dumps(dict(tensorflow=tf.__version__,
          visible_devices=[str(device) for device in tf.config.get_visible_devices()],
          intraop_threads=args.threads, interop_threads=1, input_windows=input_windows,
          requested_seq_len=settings.seq_len, effective_seq_len=effective)), flush=True)
    # Reset diagnostic state; inference state is constructed anew by upstream.
    from tiberius.main import log_config
    log_config.clear()
    # Preserve the predictor output before its own length/ORF filters.
    import bricks2marble as b2m
    raw_path = args.out + '.raw.gff3'
    with open(raw_path, 'w'):
        pass
    original_filter = b2m.tools.check_min_coding_length
    def capture_before_filter(annotation, *values, **options):
        annotation.to_gff3(raw_path, source='Tiberius', append=True)
        return original_filter(annotation, *values, **options)
    b2m.tools.check_min_coding_length = capture_before_filter
    from tiberius import eval_model_class as evaluation
    original_load = evaluation.PredictionGTF.load_model
    original_factory = evaluation.lstm_model
    timings = {}
    def load_model(instance, *values, **options):
        started = time.perf_counter()
        key = instance.model_path
        reusable = state is not None and state.get('path') == key
        if reusable:
            network = state['network']
            original_weights = network.load_weights
            evaluation.lstm_model = lambda **kwargs: network
            network.load_weights = lambda *a, **k: None
        try:
            original_load(instance, *values, **options)
            if state is not None:
                state.update(path=key, network=instance.lstm_model)
            timings.update(model_load_seconds=time.perf_counter()-started, neural_weights_reused=reusable)
        finally:
            evaluation.lstm_model = original_factory
            if reusable:
                network.load_weights = original_weights
    evaluation.PredictionGTF.load_model = load_model
    started = time.perf_counter()
    try:
        run_tiberius(settings)
    finally:
        b2m.tools.check_min_coding_length = original_filter
        evaluation.PredictionGTF.load_model = original_load
    timings["prediction_seconds"] = time.perf_counter()-started
    print("GENE_CONTEXT_TIMING="+json.dumps(timings), flush=True)


def initialize(args):
    os.environ.update(CUDA_VISIBLE_DEVICES="-1", TF_NUM_INTRAOP_THREADS=str(args.threads),
                      TF_NUM_INTEROP_THREADS="1", OMP_NUM_THREADS=str(args.threads),
                      OPENBLAS_NUM_THREADS=str(args.threads), TF_CPP_MIN_LOG_LEVEL="1")
    import tensorflow as tf
    tf.config.set_visible_devices([], "GPU")
    tf.config.threading.set_intra_op_parallelism_threads(args.threads)
    tf.config.threading.set_inter_op_parallelism_threads(1)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--serve", action="store_true")
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
    if min(args.threads, args.batch_size) < 1 or (not args.serve and (not args.model or not args.fasta or not args.out)):
        parser.error("model, fasta, out and positive thread/batch counts are required")
    if not args.serve:
        initialize(args)
        predict(args)
        return
    # A single resident neural network, bounded independently of window count.
    # HMMs, effective sequence lengths and upstream PredictionGTF are fresh.
    with redirect_stdout(sys.stderr):
        initialize(args)
    state = {}
    for line in sys.stdin:
        try:
            request = json.loads(line)
            settings = argparse.Namespace(**vars(args))
            for key in ['model', 'fasta', 'out']:
                setattr(settings, key, request[key])
            with open(request['stdout'], 'w') as out, open(request['stderr'], 'w') as err:
                with redirect_stdout(out), redirect_stderr(err):
                    predict(settings, state)
            response = dict(ok=True)
        except (Exception, SystemExit) as error:
            traceback.print_exc(file=sys.stderr)
            response = dict(ok=False, error=f'{type(error).__name__}: {error}')
            state.clear()
        print(json.dumps(response), flush=True)


if __name__ == "__main__":
    main()
