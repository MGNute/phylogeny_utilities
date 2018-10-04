#!/usr/bin/env python
# encoding: utf-8
 
import smtplib, sys
from datetime import datetime
from noticeEmail_params import params
 
def noticeEMail(usr, psw, fromaddr, toaddr, subj, mymsg):
    """
    Sends an email message through GMail once the script is completed.  
    Developed to be used with AWS so that instances can be terminated 
    once a long job is done. Only works for those with GMail accounts.
    
    starttime : a datetime() object for when to start run time clock

    usr : the GMail username, as a string

    psw : the GMail password, as a string 
    
    fromaddr : the email address the message will be from, as a string
    
    toaddr : a email address, or a list of addresses, to send the 
             message to
    """
 
    # Calculate run time
    # runtime=datetime.now() - starttime
    
    # Initialize SMTP server
    server=smtplib.SMTP('smtp.gmail.com:587')
    server.starttls()
    server.login(usr,psw)
    
    # Send email
    senddate=datetime.strftime(datetime.now(), '%Y-%m-%d')
    subject=subj
    m="Date: %s\r\nFrom: %s\r\nTo: %s\r\nSubject: %s\r\nX-Mailer: My-Mail\r\nContent-Type: text/html\r\n\r\n" % (senddate, fromaddr, toaddr, subject)
    msg='''
The following message has been sent:<br><br>
'''
    msg+=mymsg
    
    server.sendmail(fromaddr, toaddr, m+msg)
    server.quit()
 
def print_help():
    if len(sys.argv) == 1:
        o = '''
    No arguments detected, please add at least one:
        
        
        '''
    else:
        o = '''
        help option detected:
        
        '''

    o+='''
    A python script to send myself an email with contents.
    
    Usage:
    
    noticeEmail.py [-h (--help) | -si (--stdin)] <subject> <message>
    
    Logic: 
        - if '-h' or '--help' are included in the arguments, the script will
            print this help menu and then terminate.
        - otherwise, if '-si' or '--stdin' are included in the arguments, then
            the message of the email will be taken from <stdin>. 
            
            If this option is included, the one first argument that is *not* 
            '-si' or '--stdin' is used as the subject line for the email. Any
            arguments thereafter are ignored.
            
            If no argument is included other than '-si' and automatic subject
            line is generated included in the start and end datetime.
            
        - otherwise, the first argument is used as the subject and the second
            argument (if provided) is used as the message. If no second argument
            is provided, the message is blank. 
    '''
    print(o)


if __name__ == '__main__':
    nargs=len(sys.argv)
    # print('# arguments: %s' % nargs)
    # print('Usage: ' + ' '.join(sys.argv))
    # print('\n%s\n' % str(set(sys.argv).intersection(set(['-si','--stdin']))))
    # print('\n%s\n' % len(set(sys.argv).intersection(set(['-si', '--stdin']))))
    subj = ''
    mymsg = ''
    if len(sys.argv)<2 or len(set(sys.argv).intersection(set(['-h','--help'])))>0:
        #print the help menu and quit
        print_help()
        sys.exit(0)

    elif len(set(sys.argv).intersection(set(['-si','--stdin'])))>0:
        # mesage coming from stdin
        if nargs==2:
            # generate subject line from times
            subj+='Start: %s ----- End: ' % datetime.now()
        else:
            # get subject line from args
            i=1
            while i<nargs:
                if sys.argv[i] not in ['-si','--stdin']:
                    subj+=sys.argv[i]
                    break # ignore the rest
                i+=1

        # get the message from stdin
        # print('getting from stdin:\n')
        for ln in sys.stdin:
            if len(ln.strip())>0:
                # mymsg+=ln.strip() + '\n'
                mymsg += ln.strip() + '<br>'
            else:
                break

        # if we are generating the subject line, note the end time:
        if nargs==2:
            subj+=str(datetime.now())
    else:
        # subject is first arg, message is second (if it exists
        subj = sys.argv[1]
        if nargs>2:
            mymsg = sys.argv[2]

    params['subj']=subj
    params['mymsg']=mymsg
    noticeEMail(**params)
