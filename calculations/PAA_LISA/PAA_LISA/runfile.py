from imports import *
import plotfile
import save_fig
import writefile

def do_run(input_param):
    for keys in input_param.keys():
        if keys == 'dir_savefig':
            input_param[keys] = input_param[keys]+'/'
        globals()[keys] = input_param[keys]
        #print(vars()[keys])
        #print(dir_orbits)
    filename_list=[]
    filename_done=[]
    PAA_res={}
    count=0

    for (dirpath, dirnames, filenames) in os.walk(dir_orbits):
        print(filenames)
        for i in filenames:
            if i.split('.')[-1]=='txt':
                a = dirpath+'/'+i
                a = a.replace('//','/')
                filename_list.append(a)

    data_all={}

    for i in filename_list:
        filename_name = i.split('/')[-1]
        if i == filename_list[0]:
            new_folder=False # Adjust if you (don't) want to override
        else:
            new_folder=False
        print('Dir_extr:'+dir_extr)
        if select == 'all':
            if '/try/' in i:
                execute = False
            else:
                execute = True
        else:
            if select in i:
                execute = True
            else:
                execute = False

        if filename_name in filename_done:
            execute = False

        if execute == True:
            filename_save = i.split('/')[-1].split('_')[0]
            [data,PAA_res[filename_save]]=PAA(home = home,filename = i,directory_imp=False,read_max = length_calc,plot_on=True,dir_extr=dir_extr,new_folder=new_folder,timeunit=timeunit,LISA=LISA_opt,arm_influence=arm_influence,tstep=tstep,delay=delay,method=method,valorfunc='Function',dir_savefig=dir_savefig).PAA_func()
            PAA_res[str(count+1)] = PAA_res[filename_save]
            filename_done.append(filename_name)
            count=count+1

            data = plotfile.do_plot(data,dir_extr,i,new_folder,tstep,plot_on=plot_on)
            data = writefile.do_writefile(data,data_use=True)

            save_fig.do_save_fig(data)
            
            data_all[filename_save] = data
            data_all[str(count)] = data

    return data_all,PAA_res
