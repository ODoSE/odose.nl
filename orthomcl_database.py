#!/usr/bin/env python
"""Module to create, configure and dispose separate database instances for individual OrthoMCL runs."""

from ConfigParser import SafeConfigParser
from divergence import resource_filename
import MySQLdb
import logging as log
import os.path
import random

def _get_root_credentials():
    """Retrieve MySQL credentials from orthomcl.config to an account that is allowed to create new databases."""
    config = SafeConfigParser()
    config.read(resource_filename(__name__, 'orthomcl.cfg'))
    host = config.get('mysql', 'host')
    port = config.getint('mysql', 'port')
    user = config.get('mysql', 'user')
    passwd = config.get('mysql', 'pass')
    return host, port, user, passwd

def create_database():
    """Create database orthomcl_{random suffix}, grant rights to orthomcl user and return """
    dbname = 'orthomcl_' + str(random.getrandbits(16))
    host, port, user, passwd = _get_root_credentials()
    db_connection = MySQLdb.connect(host = host, port = port, user = user, passwd = passwd)
    cursor = db_connection.cursor()
    cursor.execute('CREATE DATABASE ' + dbname)
    cursor.execute('GRANT SELECT,INSERT,UPDATE,DELETE,CREATE VIEW,CREATE,INDEX,DROP on {0}.* TO orthomcl@{1};' \
                   .format(dbname, host))
    cursor.execute('set password for orthomcl@{host} = password(\'pass\');'.format(host = host))
    db_connection.commit()
    cursor.close()
    db_connection.close()
    log.info('Created database %s as %s on %s', dbname, user, host)
    return dbname

def get_configuration_file(run_dir, dbname):
    """Get OrthoMCL configuration file for generated database."""
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
evalueExponentCutoff=-5
oracleIndexTblSpc=NONE""".format(dbname = dbname, host = host, port = port)

    #Write to file & return file
    config_file = os.path.join(run_dir, 'orthomcl_{0}.cfg'.format(dbname))
    with open(config_file, mode = 'w') as write_handle:
        write_handle.write(config)
    return config_file

def delete_database(dbname):
    """Delete database after running OrthoMCL analysis."""
    host, port, user, passwd = _get_root_credentials()
    db_connection = MySQLdb.connect(host = host, port = port, user = user, passwd = passwd)
    cursor = db_connection.cursor()
    cursor.execute('DROP DATABASE ' + dbname)
    db_connection.commit()
    cursor.close()
    db_connection.close()
    log.info('Deleted database %s as %s from %s', dbname, user, host)
