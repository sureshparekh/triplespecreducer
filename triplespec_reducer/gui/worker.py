"""
Running a step off the interface thread.

A reduction step takes seconds to minutes; run on the window's own thread
it freezes the window, and a frozen window cannot say what it is doing.
Each step runs here, reports its progress as a line of text, and hands
back its result or the full traceback of what went wrong.
"""

import traceback

from PySide6.QtCore import QThread, Signal


class Job(QThread):
    progress = Signal(str)
    done = Signal(object)
    failed = Signal(str)

    def __init__(self, fn, *args, **kwargs):
        super().__init__()
        self.fn, self.args, self.kwargs = fn, args, kwargs

    def say(self, *parts):
        self.progress.emit("  ".join(str(p) for p in parts if p))

    def run(self):
        try:
            self.done.emit(self.fn(self.say, *self.args, **self.kwargs))
        except Exception as exc:          # every failure goes back to the window
            self.failed.emit(f"{exc}\n\n{traceback.format_exc()}")
