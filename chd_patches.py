# Fixes for chromatinhd 0.4.3, imported by steps 3-5 right after chromatinhd:
#   import chromatinhd as chd
#   import chd_patches
#
# torch_scatter: upstream commit 04086633 ("remove torch_scatter dependency")
# commented out `import torch_scatter`, but the pred models (models/pred/model/
# better.py and shared.py) still call torch_scatter.segment_sum_coo /
# segment_mean_coo when pooling fragment embeddings per cell x region, so
# training fails with NameError. This gives those modules a pure-torch
# replacement. Remove once chromatinhd fixes it.
#
# polyptich panels: chromatinhd's plot panels (e.g. models.pred.plot.Pileup and
# Predictivity, used by step 5) subclass polyptich.grid.Panel and draw on
# `self.ax`. In polyptich (0.0.9 - 0.0.20 checked) Panel is Ax2, which is itself
# the matplotlib Axes placed in the figure, and sets `self.ax` to an Axes that is
# never initialised -> AttributeError: 'Axes' object has no attribute 'xaxis'.
# Make `ax` refer to the panel itself (what polyptich positions and iterates).
import types

import torch
import polyptich.grid.panel
import chromatinhd.models.pred.model.better
import chromatinhd.models.pred.model.shared


def segment_sum_coo(src, index, dim_size):
    """Sum the rows of src that share an index: out[i] = src[index == i].sum(0)."""
    out = torch.zeros((dim_size,) + src.shape[1:], dtype=src.dtype, device=src.device)
    return out.index_add_(0, index.to(torch.int64), src)


def segment_mean_coo(src, index, dim_size):
    """Mean of the rows of src that share an index; 0 where an index has no rows."""
    index = index.to(torch.int64)
    total = segment_sum_coo(src, index, dim_size)
    count = torch.zeros(dim_size, dtype=src.dtype, device=src.device)
    count.index_add_(0, index, torch.ones_like(index, dtype=src.dtype))
    return total / count.clamp(min=1).view((-1,) + (1,) * (src.dim() - 1))


torch_scatter = types.SimpleNamespace(segment_sum_coo=segment_sum_coo, segment_mean_coo=segment_mean_coo)
chromatinhd.models.pred.model.better.torch_scatter = torch_scatter
chromatinhd.models.pred.model.shared.torch_scatter = torch_scatter


polyptich.grid.panel.Ax2.ax = property(lambda self: self, lambda self, value: None)
