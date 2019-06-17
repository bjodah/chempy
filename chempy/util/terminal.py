from contextlib import contextmanager
import os
import shutil
import sys
import time
import logging
import inspect
import pprint
import subprocess
import textwrap


class Echo:
    """ Context maganger for echoing variable assignments (in CPython) """
    def __init__(self, msg, indent='  '):
        self.msg = msg
        self.indent = indent
        self.parent_frame = inspect.currentframe().f_back

    def __enter__(self):
        print(self.msg)
        self.locals_on_entry = self.parent_frame.f_locals.copy()

    def __exit__(self, exc_t, exc_v, tb):
        new_locals = dict((k, v) for k, v in self.parent_frame.f_locals.items() if k not in self.locals_on_entry)
        print(textwrap.indent(pprint.pformat(new_locals), self.indent))


class Notify:
    def __init__(self, msg=lambda t: "Job finished in %.1f seconds" % t):
        self.msg = msg
        self.t0 = time.time()

    def __enter__(self):
        pass

    def notify(self, title, message):
        if os.environ.get('DISPLAY', ''):
            subprocess.call(['notify-send', title, message])
        else:
            print(title)
            print(message)

    def __exit__(self, exc_t, exc_v, tb):
        if exc_t is None:
            title = "Success"
        else:
            title = "Failure"
        self.notify(title, self.msg(time.time() - self.t0))


class c:

    ok = '\033[92m'    # green
    fail = '\033[91m'  # red
    endc = '\033[0m'   # reset


class Timed:
    """ Utility function for timing portions of your python script

    Parameters
    ----------
    msg : str
    timer : callable
        You can switch to e.g. ``time.process_time`` but note that if other
        programs are called from e.g. ``subprocess.Popen`` the time spent
        in those subprocesses will not be included.

    Examples
    --------
    >>> t = Timed("Counting stars...").tic(); stars.count(); t.toc_and_print()  # doctest: +SKIP
    Counting stars...                                          (42.0 s) [  ok]
    >>> with Timed("Counting sheep..."):  # doctest: +SKIP
    ...     n_sheep = animals.count('sheep')
    ...
    Counting sheep...                                          (17.2 s) [  ok]

    """

    counting = False

    def __init__(self, msg=None, timer=time.time, fmt_s='.1f', out=sys.stdout):
        self.msg = msg
        self.out = out
        self.fmt_s = fmt_s
        sys.stdout.flush()
        self.timer = timer

    def tic(self):
        if self.msg is not None:
            self.out.write(self.msg)
            self.out.flush()
        self.counting = True
        self.t = self.timer()
        return self  # Allows t = Timed().tic(); integrals.calc(); t = t.toc()

    def toc(self, ok=True):
        if self.counting:
            t = self.timer() - self.t
            if self.msg is not None:
                if ok:
                    status = 'ok'
                    color = c.ok
                else:
                    status = 'error'
                    color = c.fail
                self.out.write('%{}s\n'.format(shutil.get_terminal_size()[0] - len(self.msg)) %
                               ('(%{fmt_s} s) [{c}%5s{r}]'.format(fmt_s=self.fmt_s, c=color, r=c.endc) % (t, status)))
                self.out.flush()
            return t
        else:
            raise ValueError("Not counting, did you forget to call ``.tic()`` method?")

    def __enter__(self):
        self.tic()

    def __exit__(self, exc_type, exc_value, traceback):
        self.toc(not exc_type)


@contextmanager
def limit_logging(max_lvl=logging.CRITICAL):
    """ Contextmanager for silencing logging messages.

    Examples
    --------
    >>> with limit_logging():
    ...     logger.info("you won't see this...")  # doctest: +SKIP

    """
    _ori = logging.root.manager.disable
    logging.disable(max_lvl)
    try:
        yield
    finally:
        logging.disable(_ori)
