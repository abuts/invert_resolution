#!/usr/bin/python
import os
import sys
import subprocess


def script_help():
    print "Script to submit data file to script repository"
    print "Usage: "
    print ">>commit \"fileName\" \"commit message\" "

if len(sys.argv) < 2:
    script_help();

file=sys.argv[1];
if file.lower() == '-h':
    script_help()
    sys.exit();
  
if not os.path.exists(file):
    raise KeyError(" Can not find file \""+file+"\"");
if not os.path.isfile(file):
    raise KeyError(" The "+file+" is actually a folder");

if len(sys.argv)<3:
    print " Provide commit message to remember the purpose/scope of the file you want to commit"
    commit_message= raw_input();
else:
    commit_message =''
    for i in xrange(2,len(sys.argv)):
        commit_message = commit_message+sys.argv[i]+" ";

if commit_message[0:3] != "-m\"":
    commit_message = "-m\""+commit_message+"\""

print "Commit message: ",commit_message
err=subprocess.call(["svn","update"])
if err != 0:
    raise RuntimeError(" The \"svn update\" command has not been sucsesfull\n Try to run it manually")
err = subprocess.call(["svn","add",file])
if err != 0:
    raise RuntimeError(" The \"svn add "+file+" \" command has not been sucsesfull\n Try to run this command manually and analyze what went wrong")

err = subprocess.call(["svn","commit",commit_message])
if err != 0:
    raise RuntimeError(" The \"svn commit "+commit_message+" \" command has not been sucsesfull\n Try to run this command manually and analyze what went wrong")



