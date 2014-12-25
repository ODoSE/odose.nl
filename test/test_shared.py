import os
import shutil
import tempfile
import unittest

import shared


class Test(unittest.TestCase):

    def test_extract_archive_of_files(self):
        '''
        Extract a sample zipfile and assert that the contents from the first and only file match expecttations.
        '''
        target_dir = tempfile.mkdtemp()
        try:
            archive_file = shared.resource_filename(__name__, 'data/shared/sample.txt.zip')
            files = shared.extract_archive_of_files(archive_file, target_dir)
            with open(files[0]) as reader:
                contents = reader.read()
                self.assertEqual("12345", contents)
        finally:
            shutil.rmtree(target_dir)

    def test_extract_archive_of_files_empty_zipfile(self):
        '''
        Test that we raise an assertion when trying to extract files from an empty zipfile, as ofter reported via mail.
        '''
        fakefile = tempfile.mkstemp()[1]
        target_dir = tempfile.mkdtemp()
        try:
            self.assertRaises(AssertionError, shared.extract_archive_of_files, fakefile, target_dir)
        finally:
            os.remove(fakefile)
            shutil.rmtree(target_dir)
