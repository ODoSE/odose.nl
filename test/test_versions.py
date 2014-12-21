import logging
import os.path
import unittest

from versions import SOFTWARE_DIR
import versions


class Test(unittest.TestCase):

    def setUp(self):
        self.longMessage = True
        logging.root.setLevel(logging.DEBUG)

    @unittest.skipUnless(os.path.isdir(SOFTWARE_DIR), 'SOFTWARE_DIR must be available')
    def test_main(self):
        versions.main()
