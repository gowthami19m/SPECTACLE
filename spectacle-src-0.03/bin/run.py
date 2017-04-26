import dp_error_count
import sys

ref = sys.argv[1]
orig = sys.argv[2]
corrected = sys.argv[3]

nyy, yyn, nyn,nnn, yny= dp_error_count.main_func(ref,orig,corrected)

#print "yyn"
print yyn
#print "yny"
print yny
#print "nyy"
print nyy
#print "nyn"
print nyn
#print "nnn"
print nnn
 
