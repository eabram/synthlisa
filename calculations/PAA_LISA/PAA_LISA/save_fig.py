from imports import *

def do_save_fig(self):
    figs = self.figs
    titles=['-PAA','-diffPAA','-PAA_all','-Breathing_angles','-send_receive_angle','-Wobbling_angle','-Armlengths','-Velocity','-Retarded_time','-Relative_velocity']
    if self.delay==True:
        titles.append('-PPA_estimate')
    for i in range(0,len(figs)):
        title=self.filename_save+titles[i]+'.png'
        figs[i].savefig(self.dir_savefig+title)

        print('Figure '+title+' saved in:')
        print(self.dir_savefig)



    print('')
    print('')
    plt.close()
    
