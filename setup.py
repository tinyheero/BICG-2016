import os
import sys

conf_filename = 'paths.conf'

if os.path.exists(conf_filename):
    print 'Configuration file ' + conf_filename + ' exists'
    sys.exit()

with open(conf_filename, 'w') as conf_file:
    conf_file.write('# CBW Tutorial configuration file\n\n')
    conf_file.write('# Root directory for tutorial data\n')
    conf_file.write('tutorial_directory = \n')

print 'Edit configuration file ' + conf_filename + ' to complete the setup.'