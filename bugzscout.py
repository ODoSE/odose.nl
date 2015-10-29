#!/usr/bin/env python
'''
Module to submit errors to email automatically.

@author: tbeek
'''

import logging
import os
import socket
import sys
import traceback

__author__ = "Tim te Beek"
__copyright__ = "Copyright 2011, Netherlands Bioinformatics Centre"
__license__ = "MIT"


def report_error_to_email():
    '''From an except-clause retrieve details about the error and post the results to email'''
    desc, lines = get_error_trace_lines()
    send_as_email(subject=desc, content=lines)
    return desc


def get_error_trace_lines():
    '''For the current error return an ID line & string containing a complete stacktrace and local variables.'''
    err, value, trace_back = sys.exc_info()
    logging.error(err)
    lines = []

    # Add information identifying the system
    lines.extend(os.uname())
    lines.extend(str(item) for item in socket.gethostbyaddr(socket.gethostname()) if item)
    lines.append('\n')

    while trace_back:
        # Walk through trace_back to get all frames
        frame = trace_back.tb_frame
        trace_back = trace_back.tb_next

        # Identify frame
        lines.append("%s @ %s:%s" % (frame.f_code.co_name,
                                     frame.f_code.co_filename,
                                     frame.f_lineno))
        # Get locals in frame
        for key in sorted(frame.f_locals.keys()):
            # Calling str() on an unknown object could cause an error we don't want: Prevent that.
            try:
                value = str(frame.f_locals[key])
            except:
                value = "<ERROR CALLING STR() ON %s>" % key  # pylint: disable=W0702
            # Trim to keep the email readable
            if 1000 < len(value):
                value = value[:1000] + "[trimmed]"
            lines.append("\t%20s = %s" % (key, value))
        # Separator between this frame and the next
        lines.append('\n')

    # Append nicely formatted stack trace
    lines.extend(traceback.format_exc().splitlines())

    # Get error desciption as Type(values) @ file:line
    description = '{0} @ {1}:{2}'.format(err.__name__, frame.f_code.co_filename, frame.f_code.co_name)
    return description, '\n'.join(lines)


def send_as_email(subject='', content=''):
    '''
    Send an email to admin
    :param subject: email subject
    :type subject: str
    :param content: email content
    :type content: str
    '''
    from email.mime.text import MIMEText
    import smtplib
    try:
        me = 'admin@odose.nl'
        you = 'odose@bioinformatics.nl'
        msg = MIMEText(content)
        msg["From"] = me
        msg["To"] = you
        msg["Subject"] = subject
        s = smtplib.SMTP('192.168.1.14')
        s.sendmail(me, [you], msg.as_string())
        s.quit()
    except:
        logging.error('Failed to send email: %s', subject)
        # Do not mask original error by an error in sending email
        return


if __name__ == '__main__':
    try:
        avar = 12
        alist = ['a', 1, avar]
        aobject = object()
        raise Exception('bla')
    except Exception as e:
        report_error_to_email()
        raise
