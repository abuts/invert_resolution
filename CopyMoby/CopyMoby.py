""" Simple script to extract all mobi files from calibre library and copy them in separate directory to copy them on Kindle
"""

    import os
    import shutil
    source =r"d:\users\abuts\Documents\Calibre Library"
    targ = r"d:\users\abuts\Documents\mobi"

    nummobi=0
    for dirname, dirnames, filenames in os.walk(source):
        # print path to all subdirectories first.
        #for subdirname in dirnames:
        #    print(os.path.join(dirname, subdirname))

        # print path to all filenames.
        for filename in filenames:
            fbase,fext = os.path.splitext(filename)
            if fext == '.mobi':
                #print 'found ',filename
                nummobi+=1
                source_file = os.path.join(dirname,filename)
                targ_file = os.path.join(targ,filename)
                if not os.path.exists(targ_file):
                    shutil.copyfile(source_file,targ_file)
    print "Prepared {0} files".format(nummobi)

        # Advanced usage:
        # editing the 'dirnames' list will stop os.walk() from recursing into there.
        #if '.git' in dirnames:
            # don't go into any .git directories.
         #   dirnames.remove('.git')