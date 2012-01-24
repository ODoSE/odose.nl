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

URL_JOBS = 'https://ws2.grid.sara.nl/apps/prod/jobstates/'
URL_APPS = 'https://ws2.grid.sara.nl/apps/prod/applications/'

LSGP_USER = ''
LSGP_PASS = ''


def _load_credentials():
    """Read Life Science Grid Portal credentials from file and store in global variables."""
    #Get path to credential file
    from __init__ import resource_filename
    lsgp_credentials_file = resource_filename(__name__, 'credentials/lsg-portal.cfg')

    #Parse credential file
    from ConfigParser import SafeConfigParser
    parser = SafeConfigParser()
    parser.read(lsgp_credentials_file)
    defaults = parser.defaults()

    #Overwrite global values with read values
    global LSGP_USER, LSGP_PASS
    LSGP_USER = defaults['username']
    LSGP_PASS = defaults['password']

_load_credentials()


def run_application(application, params=None, files=None):
    """
    Run application with provided parameters and files.
    @param application: the application part of the Life Science Grid Portal url
    @param params: dictionary mapping keys to values for use as parameters
    @param files: dictionary mapping keys to files for use as parameters
    """
    # Invoke greeter version 1.0 with the below argument
    url_greeter = URL_APPS + application
    url_job = _send_request(url_greeter, params=params, files=files)
    jobid = url_job.split('/')[-1]
    logging.debug('Submitted %s run at %s', application, url_job)

    # Wait for job to finish
    _wait_for_job(jobid)
    logging.debug('%s finished', jobid)

    # Retrieve job result file & Extract files from .tgz
    directory = _save_job_result(jobid)
    logging.debug('%s saved to %s', jobid, directory)

    # Delete job result
    _send_request(url_job, method='DELETE')
    logging.debug('%s deleted', jobid)

    return directory


def _send_request(url, params=None, files=None, method=None):
    """
    Send a request to the SARA Life Science Grid Portal, using the provided credentials to login rather than x509.
    @param url: url to request / submit to
    @param values: dictionary of data that should be POSTed to the url (optional)
    @param method: string HTTP method (optional: POST when data is provided, GET otherwise)
    """
    #Add HTTP Basic Authentication
    password_manager = urllib2.HTTPPasswordMgr()
    password_manager.add_password('Grid Portal', 'ws2.grid.sara.nl', LSGP_USER, LSGP_PASS)
    auth_handler = urllib2.HTTPBasicAuthHandler(password_manager)
    #Create opener using our above auth_handler, and the StreamingHTTPSHandler from poster
    opener = urllib2.build_opener(auth_handler, StreamingHTTPSHandler)

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
        response = opener.open(request, timeout=10)
    except HTTPError as e:
        print e
        for k in sorted(e.hdrs.keys()):
            print k, e.hdrs[k]
        print e.fp.read()
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
        jobstates = _send_request(URL_JOBS)
        jobstates = dict(line.split('\t') for line in jobstates.strip().split('\r\n')[1:])
        if jobid not in jobstates or jobstates[jobid] != 'Queued':
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
    content = _send_request(URL_JOBS + jobid)
    tempdir = tempfile.mkdtemp('_' + jobid, 'lsgp_jobid_')
    outputfile = os.path.join(tempdir, 'out.tgz')
    with open(outputfile, mode='wb') as write_handle:
        write_handle.write(content)

    #Extract all files to tempdir
    tar = tarfile.open(outputfile)
    tar.extractall(path=tempdir)
    tar.close()
    return tempdir


if __name__ == '__main__':
    DIR = run_application('greeter/1.0', {'name': 'Tim'})
    print DIR, os.listdir(DIR)
    with open(os.path.join(DIR, os.listdir(DIR)[-1])) as read_handle:
        print read_handle.read()
