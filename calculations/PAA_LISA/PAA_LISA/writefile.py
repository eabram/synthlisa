from imports import *

def do_writefile(self,data_use=False,folder='output_files',filename='',orbit_name='',variable_dict={},mark=';   '):
    if data_use==True:
        
        folder = self.dir_savefig
        variable_dict = self.write_info
    
    
    if orbit_name=='':
        orbit_name = self.filename_save
    
    folder = folder.replace('/figures/','/output_files/')+orbit_name+'/'

    if not os.path.exists(folder):
        os.makedirs(folder)         

    if filename!='':
        filename = filename+'-'

    for key in variable_dict.keys():
        title = filename+key+'.txt'
        title = title.replace('/',' by ')
        title = folder+'/'+title
        writefile = open(title,'w')
        for j in range(0,len(variable_dict[key][0])):
            line=''
            for i in range(0,len(variable_dict[key])):
                line=line+str(variable_dict[key][i][j])+mark
            line = line+'\n'
            writefile.write(line)
        writefile.close()
        
        print('Data '+title+' saved in:')
        print(folder)
        print('')

    return self




