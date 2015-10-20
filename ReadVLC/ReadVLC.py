
from urllib2 import *
import sys
import re
import base64
from urlparse import urlparse


class test_connect(object):
   def __init__(self, cfg):
       self.cfg = cfg
   def getHost(self):
       return "192.168.1.2"
   def getHttpPort(self):
       return 8080;
   def autentification(self,err,request):
       authline = err.headers['www-authenticate']
       authobj = re.compile(
                  r'''(?:\s*www-authenticate\s*:)?\s*(\w*)\s+realm=['"]([^'"]+)['"]''',
                  re.IGNORECASE)
       # this regular expression is used to extract scheme and realm
       matchobj = authobj.match(authline)
       if not matchobj:
       # if the authline isn't matched by the regular expression
       # then something is wrong
         raise IOError,' Badly formatted authentication request'

       scheme = matchobj.group(1)
       realm = matchobj.group(2)
       if scheme.lower() != 'basic':
            raise IOError,'only BASIC authentication is currently supported'
       username = self.cfg['username'];
       password = self.cfg['password']

       base64string = base64.encodestring(
                '%s:%s' % (username, password))[:-1]
       authheader =  "Basic %s" % base64string
       request.add_header("Authorization", authheader)
       try:
            handle = urlopen(request)
       except IOError:
            raise ImportError,'Invalid login or password'
       #thepage = handle.read()
       return handle





   def xmlRequest(self, request, params):
        uri = "/requests/" + request + ".xml"
        if params is not None: uri = uri + "?" + urlencode(params).replace('+', '%20')
        location = "%s:%d" % (self.getHost(), self.getHttpPort())

        request = Request("http://" + location + uri)
        try:
            resp = urlopen(request)
        except IOError, e:
            if not hasattr(e, 'code') or e.code != 401:
                resp = None;
            else:
                resp = self.autentification(e,request)

        if resp is None:
            raise IOError, "No response from Server"
        xml = parse(resp)
        resp.close()
        return xml



theurl = 'http://192.168.1.2:8080'
# if you want to run this example you'll need to supply
# a protected page with your username and password

cfg = {'username':'','password':'a415buts'}

connector = test_connect(cfg)
req = Request(theurl)

reply = connector.xmlRequest("playlist",None);
try:
    handle = urlopen(req)
except IOError, e:
    # here we *want* to fail
    pass
else:
    # If we don't fail then the page isn't protected
    print "This page isn't protected by authentication."
    sys.exit(1)

if not hasattr(e, 'code') or e.code != 401:
    # we got an error - but not a 401 error
    print "This page isn't protected by authentication."
    print 'But we failed for another reason.'
    sys.exit(1)

authline = e.headers['www-authenticate']
# this gets the www-authenticate line from the headers
# which has the authentication scheme and realm in it


authobj = re.compile(
    r'''(?:\s*www-authenticate\s*:)?\s*(\w*)\s+realm=['"]([^'"]+)['"]''',
    re.IGNORECASE)
# this regular expression is used to extract scheme and realm
matchobj = authobj.match(authline)

if not matchobj:
    # if the authline isn't matched by the regular expression
    # then something is wrong
    print 'The authentication header is badly formed.'
    print authline
    sys.exit(1)

scheme = matchobj.group(1)
realm = matchobj.group(2)
# here we've extracted the scheme
# and the realm from the header
if scheme.lower() != 'basic':
    print 'This example only works with BASIC authentication.'
    sys.exit(1)

base64string = base64.encodestring(
                '%s:%s' % (cfg['username'], cfg['password']))[:-1]
authheader =  "Basic %s" % base64string
req.add_header("Authorization", authheader)
try:
    handle = urlopen(req)
except IOError, e:
    # here we shouldn't fail if the username/password is right
    print "It looks like the username or password is wrong."
    sys.exit(1)
thepage = handle.read()
