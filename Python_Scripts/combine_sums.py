import os
import glob
import re
# import pandas as pd
os.chdir("./Stanley_EQT/stanley_sum_files")

extension = 'sum'
all_filenames = [i for i in glob.glob('*.{}'.format(extension))]
print(all_filenames)

#combine all files in the list
with open ('stanley_catsum.sum','w') as outfile:
    for files in all_filenames:
        with open(files) as infile:
            outfile.write(infile.read())
        outfile.write("\n")

# with open('stanley_catsum.sum','r+') as file:
#     for line in file:
#          if re.search('\S', line):
#           file.write(line)





with open('stanley_catsum.sum', 'r+') as input, open('Stanley.sum', 'w') as output:
    for line in input:
        if re.search('\S', line):
          output.write(line)

