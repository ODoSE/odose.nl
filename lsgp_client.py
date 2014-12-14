#!/usr/bin/env python
'''
Module to talk to the SARA Life Science Grid Portal.

@author: tbeek
'''

from datetime import datetime
from poster.encode import multipart_encode, MultipartParam
from poster.streaminghttp import StreamingHTTPSHandler
import base64
import logging
import os
import shutil
import tarfile
import tempfile
import time
import urllib2

# Verified to work with: X-Portal-Version: 3603
HOSTNAME = 'apps.grid.sara.nl'
BASE_URL = 'https://' + HOSTNAME + '/'
URL_APPS = BASE_URL + 'applications/'
URL_DBS = BASE_URL + 'databases/'
URL_JOBS = BASE_URL + 'jobstates/'

# URL opener with support for multipart file upload forms
URLLIB2_OPENER = urllib2.build_opener(StreamingHTTPSHandler)


def _load_lsg_credentials():
    """Read Life Science Grid Portal credentials from file and store in os.environ variables."""
    # Get path to credential file
    from shared import resource_filename
    lsgp_credentials_file = resource_filename(__name__, 'credentials/lsg-portal.cfg')

    # Copy template config file to actual search path when file can not be found
    if not os.path.exists(lsgp_credentials_file):
        shutil.copy(lsgp_credentials_file + '.sample', lsgp_credentials_file)
        logging.info('Copied .sample file to %s', lsgp_credentials_file)

    logging.info('Credentials not found on path: Reading credentials from %s', lsgp_credentials_file)

    # Parse credential file
    from ConfigParser import SafeConfigParser
    parser = SafeConfigParser()
    parser.read(lsgp_credentials_file)
    os.environ['lsg_username'] = parser.defaults()['lsg_username']
    os.environ['lsg_password'] = parser.defaults()['lsg_password']


def submit_application_run(application, params, files):
    """
    Submit an application run with given parameters and files. Returns jobid of submitted run.
    @param application: the application part of the Life Science Grid Portal url
    @param params: dictionary mapping keys to values for use as parameters
    @param files: dictionary mapping keys to files for use as parameters
    """
    logging.info('Submitting %s run', application)
    url_application = URL_APPS + application
    logging.info('Parameters:\n%s', params)
    logging.info('Files:\n%s', files)
    url_job = send_request(url_application, params=params, files=files)
    logging.info('Result will be at: %s', url_job)
    jobid = url_job.split('/')[-1]
    return jobid


def retrieve_run_result(jobid, max_duration=None):
    """
    Retrieve results for jobid. Returns directory containing results.
    @param jobid: id of the job to wait for
    @param max_duration: the maximum number of seconds to wait for job completion
    """
    # Wait for job to finish
    wait_for_job(jobid, max_duration)
    logging.info('Finished waiting for job result %s', jobid)
    # Retrieve job result file & Extract files from .tgz
    directory = _save_job_result(jobid)
    logging.info('Saved job result %s to %s', jobid, directory)
    # Delete job result
    send_request(URL_JOBS + jobid, method='DELETE')
    logging.info('Deleted remote job result %s', jobid)
    return directory


def run_application(application, params=None, files=None, max_duration=None):
    """
    Run application with provided parameters and files, and retrieve result directly.
    @param application: the application part of the Life Science Grid Portal url
    @param params: dictionary mapping keys to values for use as parameters
    @param files: dictionary mapping keys to files for use as parameters
    """
    jobid = submit_application_run(application, params, files)
    directory = retrieve_run_result(jobid, max_duration)
    return directory


def upload_database(database, dbtype='formatdb', shared=False):
    """
    Upload a database file to the Life Science Grid Portal. Returns the url of the database file in the portal.
    @param database: database file
    @param dbtype: one of: FASTA, csbfa, formatdb
    @param shared: boolean to indicate whether this database should be shared with other users
    """
    # Build a unique URL using todays date
    today = datetime.today()
    # The trailing slash is key.. Time spent: ~2 hours
    version = str(today.date()) + '_' + str(today.time()).replace(':', '-') + '/'
    url_db_version = URL_DBS + 'www.odose.nl/' + version

    # Build up parameters for send_request
    params = {'type': dbtype,
              'shared': 1 if shared else 0}
    files = {'dbfile': database}

    # Upload
    content = send_request(url_db_version, params=params, files=files)
    logging.info('Uploaded %s database %s to: ' + content, dbtype, database)
    return content


def send_request(url, params=None, files=None, method=None):
    """
    Send a request to the SARA Life Science Grid Portal, using the provided arguments. Returns text content.
    @param url: url to request / submit to
    @param params: dictionary of parameters that should be POSTed to the url (optional)
    @param files: dictionary of files that should be POSTed to the url (optional)
    @param method: string HTTP method (optional: POST when params of files are provided, GET otherwise)
    """
    # Encode data
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

    # Create request
    headers.update(Accept='text/plain')
    request = urllib2.Request(url=url, headers=headers, data=data)

    # Set method, which could be DELETE
    if method:
        request.get_method = lambda: method

    # Add authentication, explicitly not using the urllib2.HTTPBasicAuthHandler, as it caused frequent failures
    if 'lsg_username' not in os.environ or 'lsg_password' not in os.environ:
        _load_lsg_credentials()
    base64string = base64.encodestring(os.environ['lsg_username'] + ':' + os.environ['lsg_password']).replace('\n', '')
    request.add_header("Authorization", "Basic %s" % base64string)

    # Send request over opener and retrieve response
    try:
        response = URLLIB2_OPENER.open(request, timeout=180)
    except urllib2.HTTPError as err:
        print url
        print err
        for key in sorted(err.hdrs.keys()):
            print key, err.hdrs[key]
        raise err

    # Retrieve
    content = response.read()
    response.close()
    return content


def wait_for_job(jobid, max_duration=None):
    """
    Wait for a given job to either leave the Queued status, or disappear from the job states page completely.
    @param jobid: id of the job to wait for
    @param max_duration: the maximum number of seconds to wait for job completion
    """
    duration = 0
    timetosleep = 30
    failures = 0
    while True:
        try:
            # It would be a shame to lose a reference to all jobs, so we allow for more errors when retrieving jobstates
            jobstates = send_request(URL_JOBS)
            failures = 0

            # Retrieve state for all jobs, and convert to dictionary for easier lookup
            jobstates = dict(line.split('\t')[:2] for line in jobstates.strip().split('\r\n')[1:])
            if jobid not in jobstates:
                logging.error('Life Science Grid Portal jobid %s not found in overview', jobid)
                break
            if jobstates[jobid] != 'Queued':
                break
        except urllib2.URLError as err:
            failures += 1
            # But after five consecutive failures we just plain give up
            if 5 <= failures:
                raise err

        # Raise an error when the maximum allotted time has arrived
        if max_duration:
            assert duration < max_duration, 'Job {1} overran allotted time: {0}{1}'.format(URL_JOBS, jobid)
            duration += timetosleep

        # If we're still here: Sleep for up to two minutes before trying again
        time.sleep(timetosleep)
        if timetosleep < 120:
            timetosleep += 10

    # Check job result
    if jobstates[jobid] == 'Error':
        raise IOError('Job {1} in error: {0}{1}'.format(URL_JOBS, jobid))

    # Return url with job result
    return URL_JOBS + jobid


def _save_job_result(jobid):
    """
    Save the job result gzipped tarfile and extract it's contents to a new temporary filename.
    @param jobid:
    """
    # Retrieve the produced gzipped tarfile and write it to a tempdir
    content = send_request(URL_JOBS + jobid)
    tempdir = tempfile.mkdtemp('_' + jobid, 'lsgp_jobid_')
    outputfile = os.path.join(tempdir, 'out.tgz')
    with open(outputfile, mode='wb') as write_handle:
        write_handle.write(content)

    # Extract all files to tempdir
    tar = tarfile.open(outputfile)
    tar.extractall(path=tempdir)
    tar.close()

    # Remove the out.tgz file we created above
    os.remove(outputfile)

    return tempdir


if __name__ == '__main__':
    DIR = run_application('greeter/1.0', {'name': 'Tim'})
    print DIR, os.listdir(DIR)
    with open(os.path.join(DIR, os.listdir(DIR)[0])) as read_handle:
        print read_handle.read()
