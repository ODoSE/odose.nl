#!/usr/bin/env python
"""Module to create, configure and dispose separate database instances for individual OrthoMCL runs."""

from ConfigParser import SafeConfigParser
import MySQLdb
import collections
from datetime import datetime
import os
import shutil

import logging as log
from shared import resource_filename


__author__ = "Tim te Beek"
__copyright__ = "Copyright 2011, Netherlands Bioinformatics Centre"
__license__ = "MIT"


Credentials = collections.namedtuple('Credentials', ['host', 'port', 'user', 'passwd'])


def _get_root_credentials():
    """Retrieve MySQL credentials from orthomcl.config to an account that is allowed to create new databases."""
    orthomcl_credentials_file = resource_filename(__name__, 'credentials/orthomcl.cfg')

    # Copy template config file to actual search path when file can not be found
    if not os.path.exists(orthomcl_credentials_file):
        shutil.copy(orthomcl_credentials_file + '.sample', orthomcl_credentials_file)
        log.info('Copied .sample file to %s', orthomcl_credentials_file)

    # Parse configuration file
    config = SafeConfigParser()
    config.read(orthomcl_credentials_file)
    host = config.get('mysql', 'host')
    port = config.getint('mysql', 'port')
    user = config.get('mysql', 'user')
    passwd = config.get('mysql', 'pass')

    # Fall back to environment value for password when available
    if passwd == 'pass' and 'mysql_password' in os.environ:
        passwd = os.environ['mysql_password']

    return Credentials(host, port, user, passwd)


def create_database():
    """Create database orthomcl_{random suffix}, grant rights to orthomcl user and return """
    # Build a unique URL using todays date
    dbname = 'orthomcl_{t.year}_{t.month}_{t.day}_at_{t.hour}_{t.minute}_{t.second}'.format(t=datetime.today())
    dbhost, port, user, passwd = _get_root_credentials()
    clhost = 'odose.nl' if dbhost not in ['127.0.0.1', 'localhost'] else dbhost
    db_connection = MySQLdb.connect(host=dbhost, port=port, user=user, passwd=passwd)
    cursor = db_connection.cursor()
    cursor.execute('CREATE DATABASE ' + dbname)
    cursor.execute('GRANT ALL on {0}.* TO orthomcl@\'{1}\' IDENTIFIED BY \'pass\';'.format(dbname, clhost))
    db_connection.commit()
    cursor.close()
    db_connection.close()
    log.info('Created database %s as %s on %s', dbname, user, dbhost)
    return dbname


def get_configuration_file(run_dir, dbname, evalue_exponent):
    """Return OrthoMCL configuration file for generated database and evalue_exponent.
    dbname - unique unused database name
    evalue_exponent - BLAST similarities with Expect value exponents greater than this value are ignored"""
    host, port = _get_root_credentials()[:2]
    config = """# OrthoMCL configuration file for generated database
dbVendor=mysql
dbConnectString=dbi:mysql:{dbname}:{host}:{port}
dbLogin=orthomcl
dbPassword=pass
similarSequencesTable=SimilarSequences
orthologTable=Ortholog
inParalogTable=InParalog
coOrthologTable=CoOrtholog
interTaxonMatchView=InterTaxonMatch
percentMatchCutoff=50
evalueExponentCutoff={evalue_exponent}
oracleIndexTblSpc=NONE""".format(dbname=dbname, host=host, port=port, evalue_exponent=evalue_exponent)

    # Write to file & return file
    config_file = os.path.join(run_dir, '{0}.cfg'.format(dbname))
    with open(config_file, mode='w') as write_handle:
        write_handle.write(config)
    return config_file


def delete_database(dbname):
    """Delete database after running OrthoMCL analysis."""
    host, port, user, passwd = _get_root_credentials()
    db_connection = MySQLdb.connect(host=host, port=port, user=user, passwd=passwd)
    cursor = db_connection.cursor()
    cursor.execute('DROP DATABASE ' + dbname)
    db_connection.commit()
    cursor.close()
    db_connection.close()
    log.info('Deleted database %s as %s from %s', dbname, user, host)
