import logging
import os
import tempfile
import unittest

import select_taxa


class Test(unittest.TestCase):

    def setUp(self):
        self.longMessage = True
        logging.root.setLevel(logging.DEBUG)

    def test_main(self):
        '''
        Select a single genome and assert the download log file contains the correct output for it.
        '''
        # Setup arguments
        target = tempfile.mktemp()[1]
        try:
            args = ('--genomes=13960 --genomes-file=' + target).split()

            # Write to argument file
            select_taxa.main(args)

            # Assert contents
            with open(target) as reader:
                contents = reader.read()
            self.assertIn('17745.1\tEscherichia coli E24377A', contents)
        finally:
            os.remove(target)
