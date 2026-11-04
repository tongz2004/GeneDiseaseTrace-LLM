"""Run RENET2's own prediction entry point under current PyTorch.

Usage: .venv/Scripts/python.exe run_predict.py <arguments of renet2.predict>

RENET2 targets PyTorch 1.2. Under current PyTorch its checkpoints need weights_only=False, and they pickle a training
scheduler whose class has changed, so this loads the model weights alone. Prediction never uses the optimizer or
scheduler; the model, weights and decision rule are RENET2's.
"""
import functools
import sys
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "vendor" / "RENET2" / "src"))
import torch  # noqa: E402

torch.load = functools.partial(torch.load, weights_only=False)
torch.optim.lr_scheduler.ReduceLROnPlateau._cmp = lambda self, *args: False  # referenced by the pickled scheduler

_pack = torch.nn.utils.rnn.pack_padded_sequence
torch.nn.utils.rnn.pack_padded_sequence = lambda x, lengths, **kw: _pack(x, lengths.cpu(), **kw)  # lengths must be on CPU now

from renet2 import predict  # noqa: E402


def load_checkpoint(config, checkpoint_file):
    model = predict.Base_Net(config)
    model.load_state_dict(torch.load(checkpoint_file, map_location="cpu")["model"])
    return model, None, None


predict.load_checkpoint = load_checkpoint

if __name__ == "__main__":
    try:
        predict.main()
    except BaseException:
        traceback.print_exc()
        raise
