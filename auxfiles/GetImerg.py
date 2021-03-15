#    #!/local/anaconda3/bin/Python3
import sys
import subprocess
import os

server = 'https://arthurhouhttps.pps.eosdis.nasa.gov'
def usage():
	print()
	print('Download imerg files for the given date')
	print()
	print('Usage: getImerg DATE')
	print(' DATE - Format is YYY-MM-DD')
	print()

def main(argv):
    
	# make sure the user provided a date
# 	if len(argv) != 2:
# 		usage()
# 		sys.exit(1)
	# make sure user gave a valid date
#	year, month, day = argv[1].split('-')
	# loop through the file list and get each file
    
    year = str('2018')
    month = str('04')
    day = str('11')
    
	#file_list = get_file_list(year,month,day)
    file_list = get_file_list('2018','04','11')
	for filename in file_list:
		get_file(filename)
        printf(filename)
		
def get_file_list(year, month, day):
	''' Get the file listing for the given year/month/day 
		using curl.
		Return list of files (could be empty).
	'''
# 	url = server + '/gpmdata/' + \
# 		'/'.join([year, month, day]) + \
# 		'/imerg/'
	url = server + '/gpmdata/' + year +'/'+month+'/'day + '/imerg/'
	cmd = 'curl -n ' + url
	args = cmd.split()
	
	process = subprocess.Popen(args,
		stdout=subprocess.PIPE,
		stderr=subprocess.PIPE)
	stdout = process.communicate()[0].decode()
		
	if stdout[0] == '<':
		print ('No imerg files for the given date')
		return []
	
	file_list = stdout.split()
	
	return file_list
		
def get_file(filename):
	''' Get the given file from arthurhouhttps using curl. '''
	url = server + filename
	cmd = 'curl -n ' + url + ' -o ' + \
		os.path.basename(filename)
	args = cmd.split()
	process = subprocess.Popen(args,		
		stdout=subprocess.PIPE,
		stderr=subprocess.PIPE)
	process.wait()  # wait so this program doesn't end
					# before getting all files

if __name__ == '__main__':
	main(sys.argv)
	
	
	
	
	
	
	
	
	
	
	
	
	
	