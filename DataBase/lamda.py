#/usr/bin/env python

from HTMLParser import HTMLParser
import IPython
import urllib

class AnchorParser(HTMLParser):
	def __init__(self):
		HTMLParser.__init__(self)
		self.loop = 0
		self.files = []
		
	def handle_starttag(self, tag, attrs):
		if tag =='a':
			for key, value in attrs:
				if key == 'href':
					if (value.find(".dat") != -1):
						print self.loop, value
						self.loop = self.loop + 1
						self.files.append(value)
						
class molpop():
	def __init__(self):
		pass
	def parse(self, file, data):
		
# File with energy levels and Einstein A coefficients
		print "Writing radiative data..."
		
		molName = file.split('.')[0]
		nLevels = int(data[5])
		molMass = float(data[3])
		
		f = open(molName+"_lamda.molecule", "w")
		f.write('Energy levels and radiative transitions of molecule\n')
		f.write(file+'\n')		
		f.write('\n')
		f.write('N. levels and molecular mass\n')
		f.write('>\n')
		f.write('%d   %f\n'%(nLevels,molMass))
		f.write('\n')
		f.write('      N        g      Energy in cm^{-1}    Level details\n')
		f.write('>\n')
		
		for i in range(nLevels):
			temp = data[7+i].split()
			f.write(temp[0]+'   '+temp[2]+'   '+temp[1]+"   '"+temp[3]+"'\n")
			
		f.write('\n')
		f.write('Einstein coefficients A_ij\n')
		f.write('Reference: LAMDA\n')
		f.write('\n')
		f.write('    i         j         A_ij in s^{-1}\n')
		f.write('>\n')
				
		
		nTransitions = int(data[7+nLevels+1])
		for i in range(nTransitions):
			temp = data[7+nLevels+3+i].split()
			f.write(temp[1]+' '+temp[2]+'  '+temp[3]+'\n')
			
		f.close()
						
		pointer = 7 + nLevels + 1 + nTransitions + 1 + 1 + 1

# File with collisional rates
		nCollPartners = int(data[pointer])
		pointer += 1
		
		print "Writing collisional data..."
		for i in range(nCollPartners):
			pointer += 1
			whichCollision = data[pointer].split()
			temp = whichCollision[1].split('-')
			
			fileName = file+'_'+temp[1]+"_lamda.kij"
			
			print "Collisions with "+temp[1]+" -> "+'Coll/'+fileName
			
			f = open('Coll/'+fileName, "w")
			
			pointer += 2
			nCollisions = int(data[pointer])
			
			pointer += 2
			nTemperatures = int(data[pointer])
			
			pointer += 2
			Temperatures = data[pointer]
			
			pointer += 2
		
			f.write('Collision rate coefficients \n')
			f.write('Reference: LAMDA\n')
			f.write('\n')
			f.write('>\n')
			f.write('\n')
			f.write('Number of temperature columns = '+str(nTemperatures)+'\n')
			f.write('\n')
			f.write('  I   J                        TEMPERATURE (K)\n')
			f.write('\n')
			f.write(Temperatures)
			f.write('\n')
			
			for j in range(nCollisions):
				f.write(data[pointer])
				pointer += 1
				
			f.close()
			

# Parse the directory with all the molecules
parser = AnchorParser()
data = urllib.urlopen('http://home.strw.leidenuniv.nl/~moldata/datafiles/').read()
print "List of available molecules"
parser.feed(data)

# Select the desired molecule
nb = raw_input("Select which one to download (separated with spaces if many) (0-"+str(parser.loop-1)+") ")

nb = nb.split()

for iterator in nb:
	indexMolecule = int(iterator)
	print "Downloading "+parser.files[indexMolecule]

# Download the molecule and parse it, generating the appropriate files
	ur = urllib.urlopen('http://home.strw.leidenuniv.nl/~moldata/datafiles/'+parser.files[indexMolecule])
	data = ur.readlines()

	mol = molpop()

	mol.parse(parser.files[indexMolecule],data)