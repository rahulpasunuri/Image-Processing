import os
import random

print "make"
print

dirname = "attraction_images"
files = os.listdir(dirname)

program = "./p4 part3.2"
allFiles=""


att = {}

for f in files:
	index = f.find('_')
	name=f[:index]
	if name in att:
		att[name].append(f)
	else:
		att[name] = [f]
	allFiles = allFiles + " "+dirname+"/"+ f

for key in att:
	#generate query for this key.
	#randomize the list..
	random.shuffle(att[key])	
	queryImage = att[key][0]
	
	print ""+program+" "+dirname+"/"+queryImage+" "+allFiles
	print
