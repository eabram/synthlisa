filename='NGO_1M_10deg_synthlisa.txt'

file_read=open(filename,'r')
file_write=open(filename.split('.')[0]+'_adjust.txt','w')

for line in file_read.readlines():
    a=line.split(' ')
    if line[0]!= '#':
        t=a[0]
        a_new=a[1:len(a)]
        a_new.append(t)
        linew=line
    else:
        a_new = a
        linew=''
        for i in a_new:
            linew=linew+i+'    '
    file_write.write(linew)
file_read.close()
file_write.close()







