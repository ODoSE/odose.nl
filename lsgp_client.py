#!/usr/bin/env python
'''
Module to talk to the SARA Life Science Grid Portal.

@author: tbeek
'''

from poster.encode import multipart_encode, MultipartParam
from poster.streaminghttp import StreamingHTTPSHandler
from urllib2 import HTTPError
import logging
import os
import tarfile
import tempfile
import time
import urllib2

URL_APPS = 'https://ws2.grid.sara.nl/apps/prod/applications/'
URL_DBS = 'https://ws2.grid.sara.nl/apps/prod/databases/'
URL_JOBS = 'https://ws2.grid.sara.nl/apps/prod/jobstates/'

URLLIB2_OPENER = None


def _build_authenticated_multipart_opener():
    """Read Life Science Grid Portal credentials from file and store in global variables."""
    #Get path to credential file
    from __init__ import resource_filename
    lsgp_credentials_file = resource_filename(__name__, 'credentials/lsg-portal.cfg')

    #Parse credential file
    from ConfigParser import SafeConfigParser
    parser = SafeConfigParser()
    parser.read(lsgp_credentials_file)
    defaults = parser.defaults()

    #Add HTTP Basic Authentication
    password_manager = urllib2.HTTPPasswordMgr()
    password_manager.add_password('Grid Portal', 'ws2.grid.sara.nl', defaults['username'], defaults['password'])
    auth_handler = urllib2.HTTPBasicAuthHandler(password_manager)

    #Create opener using our above auth_handler, and the StreamingHTTPSHandler from poster to handle multipart forms
    global URLLIB2_OPENER
    URLLIB2_OPENER = urllib2.build_opener(auth_handler, StreamingHTTPSHandler)

_build_authenticated_multipart_opener()


def submit_application_run(application, params, files):
    """
    Submit an application run with given parameters and files. Returns jobid of submitted run.
    @param application: the application part of the Life Science Grid Portal url
    @param params: dictionary mapping keys to values for use as parameters
    @param files: dictionary mapping keys to files for use as parameters
    """
    url_application = URL_APPS + application
    url_job = send_request(url_application, params=params, files=files)
    logging.info('Submitted %s run; job result will be at: %s', application, url_job)
    jobid = url_job.split('/')[-1]
    return jobid


def retrieve_run_result(jobid):
    """
    Retrieve results for jobid. Returns directory containing results.
    @param jobid:
    """
    # Wait for job to finish
    _wait_for_job(jobid)
    logging.info('Finished waiting for job result %s', jobid)
    # Retrieve job result file & Extract files from .tgz
    directory = _save_job_result(jobid)
    logging.info('Saved job result %s to %s', jobid, directory)
    # Delete job result
    send_request(URL_JOBS + jobid, method='DELETE')
    logging.info('Deleted remote job result %s', jobid)
    return directory


def run_application(application, params=None, files=None):
    """
    Run application with provided parameters and files, and retrieve result directly.
    @param application: the application part of the Life Science Grid Portal url
    @param params: dictionary mapping keys to values for use as parameters
    @param files: dictionary mapping keys to files for use as parameters
    """
    jobid = submit_application_run(application, params, files)
    directory = retrieve_run_result(jobid)
    return directory


def upload_database(database, dbtype='formatdb', shared=False):
    """
    Upload a database file to the Life Science Grid Portal. Returns the url of the database file in the portal.
    @param database: database file
    @param dbtype: one of: FASTA, csbfa, formatdb
    @param shared: boolean to indicate whether this database should be shared with other users
    """
    #Build a unique URL using todays date
    from datetime import datetime
    today = datetime.today()
    # The trailing slash is key.. Time spent: ~2 hours
    version = str(today.date()) + '_' + str(today.time()).replace(':', '-') + '/'
    url_db_version = URL_DBS + 'www.odose.nl/' + version

    #Build up parameters for send_request
    params = {'type': dbtype,
              'shared': 1 if shared else 0}
    files = {'dbfile': database}

    #Upload
    content = send_request(url_db_version, params=params, files=files)
    logging.info('Uploaded %s database %s to: ' + content, dbtype, database)
    return content


def send_request(url, params=None, files=None, method=None):
    """
    Send a request to the SARA Life Science Grid Portal, using the provided arguments. Returns text content.
    @param url: url to request / submit to
    @param values: dictionary of data that should be POSTed to the url (optional)
    @param method: string HTTP method (optional: POST when data is provided, GET otherwise)
    """
    #Encode data
    data = None
    headers = {}
    multipart_params = []
    if params:
        for key, value in params.iteritems():
            multipart_params.append(MultipartParam(key, value))
    if files:
        for key, value in files.iteritems():
            multipart_params.append(MultipartParam.from_file(key, value))
    if multipart_params:
        data, headers = multipart_encode(multipart_params)

    #Create request
    headers.update(Accept='text/plain')
    request = urllib2.Request(url=url, headers=headers, data=data)
    if method:
        request.get_method = lambda: method

    #Send request over opener and retrieve response
    try:
        response = URLLIB2_OPENER.open(request, timeout=60)
    except HTTPError as e:
        print e
        for key in sorted(e.hdrs.keys()):
            print key, e.hdrs[key]
        raise e

    #Retrieve
    content = response.read()
    response.close()
    return content


def _wait_for_job(jobid):
    """
    Wait for a given job to either leave the Queued status, or disappear from the job states page completely.
    @param jobid: id of the job to wait for
    """
    duration = 10
    while True:
        #Retrieve state for all jobs, and convert to dictionary for easier lookup
        jobstates = send_request(URL_JOBS)
        jobstates = dict(line.split('\t') for line in jobstates.strip().split('\r\n')[1:])
        if jobid not in jobstates:
            logging.error('Life Science Grid Portal jobid %s not found in overview', jobid)
            break
        if jobstates[jobid] != 'Queued':
            break
        #Sleep for up to two minutes
        time.sleep(duration)
        if duration < 120:
            duration += 10


def _save_job_result(jobid):
    """
    Save the job result gzipped tarfile and extract it's contents to a new temporary filename.
    @param jobid:
    """
    #Retrieve the produced gzipped tarfile and write it to a tempdir
    content = send_request(URL_JOBS + jobid)
    tempdir = tempfile.mkdtemp('_' + jobid, 'lsgp_jobid_')
    outputfile = os.path.join(tempdir, 'out.tgz')
    with open(outputfile, mode='wb') as write_handle:
        write_handle.write(content)

    #Extract all files to tempdir
    tar = tarfile.open(outputfile)
    tar.extractall(path=tempdir)
    tar.close()

    #Remove the out.tgz file we created above
    os.remove(outputfile)

    return tempdir


if __name__ == '__main__':
    DIR = run_application('greeter/1.0', {'name': 'Tim'})
    print DIR, os.listdir(DIR)
    with open(os.path.join(DIR, os.listdir(DIR)[-1])) as read_handle:
        print read_handle.read()
