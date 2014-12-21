import MySQLdb
import shutil
import tempfile
import unittest

import orthomcl_database


class Test(unittest.TestCase):

    def setUp(self):
        self.run_dir = tempfile.mkdtemp()
        self.credentials = orthomcl_database._get_root_credentials()

    def tearDown(self):
        shutil.rmtree(self.run_dir)

    def test_get_configuration_file(self):
        '''
        Create a configuration file, and ensure the contents match assumptions.
        '''
        conffile = orthomcl_database.get_configuration_file(self.run_dir, 'test_dbname', 5)
        with open(conffile) as reader:
            content = reader.read()
        self.assertIn('orthomcl', content)
        self.assertIn('127.0.0.1', content)
        self.assertIn('mysql', content)
        self.assertIn('evalueExponentCutoff=5\n', content)

    def test_create_database(self):
        '''
        Create a database, connect to it and perform a simple select query, verify the outcome and delete the database.
        '''
        try:
            # Create database
            dbname = orthomcl_database.create_database()

            # Access database as restricted user
            db_connection = MySQLdb.connect(host=self.credentials.host,
                                            port=self.credentials.port,
                                            user='orthomcl', passwd='pass')
            db_connection.query('SELECT 1')
            result = db_connection.store_result()
            self.assertEqual(1L, result.fetch_row()[0][0])
            db_connection.close()
        finally:
            if dbname:
                # Delete database
                orthomcl_database.delete_database(dbname)
