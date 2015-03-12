import urllib2
import sys
import magic
import os


def downLoadfile(url,savedir,counter):
	try:
		data = urllib2.urlopen(url,timeout=4)
		with open("try","wb") as f:
			f.write(data.read())
		mime = magic.Magic(mime=True)
		mimetype=mime.from_file("try")
		print mimetype
		if "image" not in mimetype:
			print "image not in mimetype"
			return False
		else:
			extension=mimetype.split("/")[1]
			saveloc=savedir+"/"+str(counter)+"."+extension
			with open(saveloc,"wb") as f:
				f.write(open("try","rb").read())
			print "saved as",saveloc
			return True
	except: # (urllib2.HTTPError,urllib2.URLError,IOError):	
		print "something wrong happened"
		return False

def main():	
	content=open(sys.argv[1]).read().split("\n")[1203:]
	counter=769
	for index,item in enumerate(content):
		print index,item
		directory=item.split("\t")[0]
		url=item.split("\t")[1]
		if not os.path.exists(directory):
    			os.makedirs(directory)
		if downLoadfile(url,directory,counter):
			counter+=1
	

if __name__=="__main__":
	main()
