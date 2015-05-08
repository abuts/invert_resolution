#!/usr/bin/python 
from xml.dom import minidom 
import smtplib
from email.mime.text import MIMEText

def checker() :
    domObj=minidom.parse(r"/opt/Mantid/instrument/LET_Parameters.xml")
    #domObj=minidom.parse(r"d:/Data/Mantid_GIT/Code/Mantid/instrument/LET_Parameters.xml")
    params = domObj.getElementsByTagName("parameter")

    val = process_param(params,"vanadium-mass");

    print "found value: ",val," of type : ",type(val)

    mm = float(val)
    print " mass = ",mm

    if val != "20.79" :
        me = "Alex.Buts@stfc.ac.uk"
        you = "Alex.Buts@stfc.ac.uk"

        msg = MIMEText("vanadium changed ")
        # me == the sender's email address
        # you == the recipient's email address
        msg['Subject'] = 'The vanadium mass is: '+val
        msg['From'] = me
        msg['To']   = you
        # Send the message via our own SMTP server, but don't include the
        # envelope header.
        s = smtplib.SMTP('localhost')
        s.sendmail(me, [you], msg.as_string())
        s.quit()
    else:
        print " Sleep tight"





def process_param(params,parName):
    for param in params :
        par_name=param.getAttribute("name")
        if par_name == parName :
             values= param.getElementsByTagName('value') 
             return values[0].attributes['val'].value 







if __name__ == '__main__' :
    print " start checker"
    checker()



