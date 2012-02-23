"""Utilities to facilitate the writing of tests for pypsy.
"""
#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

# Third-party

#-----------------------------------------------------------------------------
# Functions and classes
#-----------------------------------------------------------------------------

def test(doctests=True, first_package_wins=True, extra_argv=None):
    """

    Run the pypsy test suite using nose.

    Parameters
    ----------

    doctests: bool, optional
       Whether to run the doctests. Defaults to True

    extra_argv: string, list or tuple, optional
       Additional argument (string) or arguments (list or tuple of strings) to
       be passed to nose when running the tests.

    """
    from numpy.testing import noseclasses
    # We construct our own argv manually, so we must set argv[0] ourselves
    argv = ['nosetests',
            # Name the package to actually test, in this case pypsy
            'pypsy',

            # extra info in tracebacks
            '--detailed-errors',

            # We add --exe because of setuptools' imbecility (it blindly does
            # chmod +x on ALL files).  Nose does the right thing and it tries
            # to avoid executables, setuptools unfortunately forces our hand
            # here.  This has been discussed on the distutils list and the
            # setuptools devs refuse to fix this problem!
            '--exe',
            ]

    # If someone wants to add some other argv
    if extra_argv is not None:
        if isinstance(extra_argv, list) or isinstance(extra_argv, list):
            for this in extra_argv: argv.append(this)
        else:
            argv.append(extra_argv)

    if doctests:
        argv.append('--with-doctest')

    # Now nose can run
    return noseclasses.NumpyTestProgram(argv=argv, exit=False).result

# Tell nose that the test() function itself isn't a test, otherwise we get a
# recursive loop inside nose.
test.__test__ = False
