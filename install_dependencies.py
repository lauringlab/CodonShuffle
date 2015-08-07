import subprocess
subprocess.call(['pip', 'install', 'numpy'])
subprocess.call(['pip', 'install', 'scipy'])
subprocess.call(['pip', 'install', 'matplotlib'])
subprocess.call(['pip', 'install', 'pandas'])
subprocess.call(['pip', 'install', 'statsmodels'])
subprocess.call(['pip', 'install', 'patsy'])
subprocess.call(['pip', 'install', 'biopython'])
#subprocess.call(['easy_install', '-f', 'http://biopython.org/DIST/', 'biopython'])
subprocess.call(['pip', 'install', 'ggplot'])
print "Ready to run CodonShuffle"