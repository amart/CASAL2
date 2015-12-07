import os
import sys
import subprocess
import os.path
import shutil
import fileinput
import re
from datetime import datetime, date
from dateutil import tz

EX_OK = getattr(os, "EX_OK", 0)

class DebBuilder:
  def start(self):
    print '-- Starting Deb Builder'
    p = subprocess.Popen(['git', '--no-pager', 'log', '-n', '1', '--pretty=format:%H%n%h%n%ci' ], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    lines = out.split('\n')          
    if len(lines) != 3:
      return Globals.PrintError('Format printed by GIT did not meet expectations. Expected 3 lines but got ' + str(len(lines)))
    print '-- CASAL2 Revision: ' + lines[1]
    if not os.path.exists('bin/linux/deb'):
      os.mkdir('bin/linux/deb')

    folder = 'bin/linux/deb/CASAL2_' + lines[1]    
    os.system('rm -rf ' + folder)
    os.makedirs(folder + '/usr/local/bin')
    os.system('cp bin/linux/release/casal2 ' + folder + '/usr/local/bin')
    os.makedirs(folder + '/DEBIAN')
    control_file = open(folder + '/DEBIAN/control', 'w')
    control_file.write('Package: CASAL2\n')
    control_file.write('Version: ' + lines[1] + '\n')
    control_file.write('Section: base\n')
    control_file.write('Priority: optional\n')
    control_file.write('Architecture: amd64\n')
    control_file.write('Maintainer: Scott Rasmussen (Zaita) <scott@zaita.com>\n')
    control_file.write('Description: CASAL2 Modeling Platform\n')
    control_file.close()

    if os.system('dpkg-deb --build ' + folder) != EX_OK:
      return Globals.PrintError('Failed to build DEB Package')
    return True