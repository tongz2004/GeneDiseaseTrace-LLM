"""BioREx prediction without BioREx's TF trainer wrapper, which crashes (segmentation fault) on native Windows.

Usage: .venv/Scripts/python.exe experiments/biorex/predict.py IN_TSV OUT_DIR [--batch 16]

Reproduces the prediction path of vendor/BioREx/src/run_ncbi_rel_exp.py (--do_predict): text from column 7 of the
converter's TSV (no header), the tokenizer saved with the pretrained model (its added tokens already include the
entity tags), truncation at 512 tokens, the TF model's logits written one row per instance to OUT_DIR/test_results.tsv
in BioREx's format, so src/utils/run_pubtator_eval.py --exp_option to_pubtator can read them unchanged. Padding is to
the longest instance in the batch instead of always 512; with the attention mask this gives the same logits.
"""
import argparse
import os
import sys
from pathlib import Path

os.environ.setdefault("TF_USE_LEGACY_KERAS", "1")
os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "2")

import tensorflow  # noqa: E402,F401  (before pandas: the other order crashes TF on native Windows)
import pandas as pd  # noqa: E402
from transformers import AutoTokenizer, TFAutoModelForSequenceClassification  # noqa: E402

MODEL = Path(__file__).resolve().parent / "vendor" / "BioREx" / "pretrained_model_biolinkbert"
TEXT_COLUMN = 7


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("in_tsv")
    parser.add_argument("out_dir")
    parser.add_argument("--batch", type=int, default=16)
    args = parser.parse_args()
    texts = pd.read_csv(args.in_tsv, sep="\t", header=None, dtype=str).fillna("None").iloc[:, TEXT_COLUMN].tolist()
    tokenizer = AutoTokenizer.from_pretrained(MODEL)
    model = TFAutoModelForSequenceClassification.from_pretrained(MODEL)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    with (out / "test_results.tsv").open("w") as writer:
        for i in range(0, len(texts), args.batch):
            enc = tokenizer(texts[i:i + args.batch], truncation=True, max_length=512, padding=True, return_tensors="tf")
            for row in model(enc, training=False).logits.numpy():
                writer.write("\t".join(map(str, row)) + "\n")
            print(f"{min(i + args.batch, len(texts))}/{len(texts)}", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
