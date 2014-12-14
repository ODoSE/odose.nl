import shutil
import tempfile
import unittest

import orthomcl_database


class Test(unittest.TestCase):

    def setUp(self):
        self.run_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.run_dir)

    def test_get_configuration_file(self):
        conffile = orthomcl_database.get_configuration_file(self.run_dir, 'test_dbname', 5)
        with open(conffile) as reader:
            content = reader.read()
        self.assertIn('orthomcl', content)
        self.assertIn('127.0.0.1', content)
        self.assertIn('mysql', content)
        self.assertIn('evalueExponentCutoff=5\n', content)

    def test_create_database(self):
        dbname = orthomcl_database.create_database()
        orthomcl_database.delete_database(dbname)
