# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

try:
    import numpy as np
except ImportError:
    np = None


def broadcast_stack(*args, **kwargs):
    as_scalars = kwargs.pop('as_scalars', False)
    if kwargs != {}:
        raise ValueError("Got unknown kwargs: %s" % kwargs)
    args = [np.atleast_1d(arg) for arg in args]
    if as_scalars:
        args = [arg.reshape(arg.shape + (1,)) if arg.size > 1 else arg for arg in args]
    if all([arg.ndim == 1 for arg in args]):
        return np.concatenate(args)
    head_shape = ()
    leading_length = 0
    for arg in args:
        leading_length += arg.shape[-1]
        if arg.ndim > 1:
            if head_shape is ():
                head_shape = arg.shape[:-1]
            else:
                if arg.shape[:-1] != head_shape:
                    raise ValueError("Incompatible shapes")
    out = np.empty(head_shape + (leading_length,))
    for idx, arg in enumerate(args):
        if arg.shape[-1] != 1:
            raise ValueError("Trailing dimensions needs to be 1")
        out[..., idx] = arg[..., 0]
    return out
