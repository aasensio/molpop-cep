#/usr/bin/env python

from HTMLParser import HTMLParser
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
		f = open(file+".molecule_temp", "w")
		
		
		

# Parse the directory with all the molecules
parser = AnchorParser()
data = urllib.urlopen('http://home.strw.leidenuniv.nl/~moldata/datafiles/').read()
print "List of available molecules"
parser.feed(data)

# Select the desired molecule
nb = int(raw_input("Select which one to download (0-"+str(parser.loop-1)+") "))

print "Downloading "+parser.files[nb]

# Download the molecule and parse it, generating the appropriate files
ur = urllib.urlopen('http://home.strw.leidenuniv.nl/~moldata/datafiles/'+parser.files[nb])
data = ur.readlines()

mol = molpop()

mol.parse(parser.files[nb],data)