"""Low-level I/O helpers for silencing native libraries.

Some compiled extensions (notably ``tblite``) write directly to the process
file descriptors, so redirecting only :data:`sys.stdout` is insufficient.
The :func:`suppress_native_output` context manager performs an ``os.dup2``
swap that catches both Python and native writes for the duration of the
``with`` block.
"""

from __future__ import annotations

import contextlib
import os
import sys
from typing import Iterator


@contextlib.contextmanager
def suppress_native_output() -> Iterator[None]:
    """Redirect stdout *and* stderr (Python + native) to ``/dev/null``."""
    with open(os.devnull, "w", encoding="utf-8") as devnull:
        stdout_fd = sys.stdout.fileno()
        stderr_fd = sys.stderr.fileno()
        saved_stdout = os.dup(stdout_fd)
        saved_stderr = os.dup(stderr_fd)
        try:
            os.dup2(devnull.fileno(), stdout_fd)
            os.dup2(devnull.fileno(), stderr_fd)
            yield
        finally:
            os.dup2(saved_stdout, stdout_fd)
            os.dup2(saved_stderr, stderr_fd)
            os.close(saved_stdout)
            os.close(saved_stderr)
