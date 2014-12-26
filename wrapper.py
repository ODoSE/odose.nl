#!/usr/bin/env python
"""Module to generically wrap modules for usage from within Galaxy.
First command line argument is a module, e.g: "translate".
Second command line argument is a file to log any output to.
Further arguments are passed to the named module's main method."""

__author__ = "Tim te Beek"
__contact__ = "brs@nbic.nl"
__copyright__ = "Copyright 2011, Netherlands Bioinformatics Centre"
__license__ = "MIT"

# First argument contains fully qualified name of module to be imported
import logging
import sys
NAME = sys.argv[1]
try:
    MODULE = __import__(NAME)
except ImportError as ie:
    print('Could not import {0}'.format(NAME))
    raise

# Second argument contains name of logging output file to use
FILE_HANDLER = logging.FileHandler(sys.argv[2], mode='w')
FILE_HANDLER.setFormatter(logging.Formatter())
FILE_HANDLER.setLevel(logging.INFO)
logging.root.addHandler(FILE_HANDLER)

try:
    # Run main method within module with remaining arguments
    if sys.argv[3:]:
        MODULE.main(sys.argv[3:])
    else:
        MODULE.main()
except SystemExit:
    # Do not report SystemExit errors to FogBugz: Just exit
    raise
except AssertionError:
    # Do not report AssertionErrors to FogBugz: Not a bug we care about
    logging.exception('An assumption failed')
    raise
except:
    # Should any other error occur, report it to FogBugz automatically
    from bugzscout import report_error_to_email
    MESSAGE = report_error_to_email()
    logging.info('Automatic bug submission reported: %s', MESSAGE)
    logging.exception('An error occurred')
    raise
finally:
    # Always remove logging handler from root
    logging.root.removeHandler(FILE_HANDLER)

# Snippet to log available environment variables from inside a Galaxy tool:
# for key in sorted($searchList[2].keys())
# silent   sys.stderr.write("\t{0} = {1} ({2})\n".format(str(key), str($searchList[2][key]), type($searchList[2][key])))
# end for
