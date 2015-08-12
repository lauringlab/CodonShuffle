import subprocess
import os
subprocess.call(['pip', 'install', 'numpy']) # Install Numpy
subprocess.call(['pip', 'install', 'scipy']) # Install Scipy
subprocess.call(['pip', 'install', 'matplotlib']) # Install Matplotlib
subprocess.call(['pip', 'install', 'pandas']) # Install Pandas
subprocess.call(['pip', 'install', 'statsmodels']) # Install Statsmodels
subprocess.call(['pip', 'install', 'patsy']) # Install Patsy
subprocess.call(['pip', 'install', 'biopython']) # Install Biopython
#subprocess.call(['easy_install', '-f', 'http://biopython.org/DIST/', 'biopython'])
subprocess.call(['pip', 'install', 'ggplot']) # Install Ggplot
os.chdir('./lib/EMBOSS-6.6.0/') #Install Emboss
subprocess.call(['./configure', '--without-x']) # Install Emboss
subprocess.call(['make']) #Install Emboss
print "Ready to run CodonShuffle"

