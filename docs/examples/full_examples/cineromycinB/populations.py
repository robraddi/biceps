#  Script to be used with Chimera
usage='''

chimera --nogui --script "align.py [arg] [arg]"\n

**NOTE**
1) Loads pdb or gro into chimera
2)
3) Saves new pdb of mutated residues

'''
print(usage)


# Libraries:{{{
import os, sys, glob, re
import numpy as np
import chimera
from chimera import runCommand as rc
from Midas.midas_text import makeCommand as mc
import Midas as m
from chimera import openModels, Molecule
from chimera import dihedral
from matplotlib import pyplot as plt
#:}}}

# Methods:{{{

def load_pdb(pdb, verbose=False):
    """Loads a structure file into chimera."""

    print('Loading Structure: %s...\n'%pdb)
    mc('open %s'%pdb)

    om = chimera.openModels
    models = om.list()

    if verbose:
        # Get the information from the pdb file loaded:
        print("###############################################################")
        pdb_info = ['TITLE']#,'REMARK']
        for m in openModels.list(modelTypes=[Molecule]):
            for header in pdb_info:
                try:
                    for line in m.pdbHeaders[header]:
                        print('Model # = %s;\n Name = %s;\n%s;\n\n'%(
                            str(m.id), str(m.name), str(line)))
                except KeyError:
                    print("Model # = %s;\n Name = %s;\n Has no header... \n\n"%(
                        str(m.id), str(m.name)))
        print("There are %s models."%len(models))
        print("###############################################################\n\n")
    return models

def get_files(path):
    convert = lambda txt: int(txt) if txt.isdigit() else txt
    return sorted(glob.glob(path), key=lambda x:[convert(s) for s in re.split("([0-9]+)",x)])

def get_pops(file, nlambda):
    pops = np.loadtxt(file)[:,nlambda-1]
    return pops

def set_view():
    rc('windowsize 600 450')
    rc('background solid white')
    #rc("preset apply publication 1")
    #rc('set silhouette')
    #rc('set silhouetteColor none')
    #rc('set silhouetteWidth 1.9')
    #rc('set shadows')

def save_image(filename="img.png", dimensions=(600,450), supersample=4):
    w,h = dimensions
    units = "points"
    rc('copy file %s width %s height %s units %s \
            supersample %i '%(filename,w,h,units, supersample))

def get_dihedral(atoms):
    a1,a2,a3,a4 = atoms
    #dihedral = rc("angle #%s:@%s #%s:@%s #%s:@%s #%s:@%s"%(model, a1, model, a2, model, a3, model, a4))
    return dihedral(a1.coord(), a2.coord(), a3.coord(), a4.coord())

def quick_plot(x, y, xlabel='x', ylabel='y', name=None, Type='scatter', labels=None,
        fig_size=(12,10), xlim=None, ylim=None, show=False):

    fig = plt.figure(figsize=fig_size)
    ax = fig.add_subplot(111)

    if Type=='scatter':
        ax.scatter(x,y)

    if Type=='line':
        ax.plot(x,y,'k')
    #ax.plot(x,p(x),"k--",label="_nolegend_")

    if Type == "hist2d":
        ax.hist2d(x,y, bins=100)

    if Type == "bar":
        ax.bar(x,y)

    if labels:
        for i in range(len(labels)):
            ax.text(x[i], y[i], str(labels[i]), color='g')

    ax.set_xlabel('%s'%xlabel, fontsize=16)
    ax.set_ylabel('%s'%ylabel, fontsize=16)
    if xlim:
       ax.set_xlim(xlim[0],xlim[1])
    if ylim:
       ax.set_ylim(ylim[0],ylim[1])

    # Setting the ticks and tick marks
    ticks = [ax.xaxis.get_minor_ticks(),
             ax.xaxis.get_major_ticks()]
    marks = [ax.get_xticklabels(),
            ax.get_yticklabels()]
    for k in range(0,len(ticks)):
        for tick in ticks[k]:
            tick.label.set_fontsize(16)
    for k in range(0,len(marks)):
        for mark in marks[k]:
            mark.set_size(fontsize=16)
            mark.set_rotation(s=0)
    #fig.tight_layout()
    if name==None:
        pass
    else:
        fig.savefig('%s'%name)
    if show:
        fig.show()




#:}}}

# Main:{{{

# Get all structure files from directory
path = str(sys.argv[1])
files = get_files(path)
# Load structures in order to correspond to model numbers
[load_pdb(files[i]) for i in range(len(files))]

# Get list of models
om = chimera.openModels
models = om.list()
# align structures to reference model
for i in range(1, len(models)):
    rc("match #%s #0"%(i))

# Styles
set_view()                          # Call on your prefered rendering styles
rc("rep wire;linewidth 3.0")        # Set the representation to wire
rc("sel :; color byatom sel; ~sel") # color by heteroatom


# Get trajectories
path = str(sys.argv[2])
files = get_files(path)
nlambda = len(files)
print("Number of lambda values: %s"%nlambda)
# Get populations from file
populations_file = os.path.join(os.path.join(path.split("/")[0:-1])[0], "populations.dat")
populations = get_pops(populations_file, nlambda)
outdir = os.path.join(path.split("/")[0:-1])[0]

# Adjust transparency to match populations
for i,pop in enumerate(populations):
    pop = pop*100
    if pop < 5: # Only show if population is greater than 1%
        rc("~show #%s"%(i))
    else:
        rc("transparency %f #%s"%(float(100-pop), i))
        # set the attribute of a surface
        rc("setattr a bfactor %.2f #%s"%(pop, i))

save_image(filename="%s/top_populated.png"%outdir, dimensions=(600,450), supersample=4)

rc("~transparency") # turn transparency off
rc("rep stick")
rc("tile")
rc("labelopt info %(molecule)s | %(bfactor).2f %%")
rc("la offset 0,-7,0 @/serialNumber=1")
save_image(filename="%s/tiled_structures.png"%outdir, dimensions=(600,450), supersample=4)

# Plot the dihedral angles of our structures
dihedrals = []
labels = []
for i,pop in enumerate(populations):
    pop = pop*100
    if pop > 5:
        #dihedrals.append((get_dihedral(i,"C5","C6","C7","C8"),get_dihedral(i,"C6","C7","C8","C9")))
        rc("sel #%s:@C5 #%s:@C6 #%s:@C7 #%s:@C8"%(i,i,i,i))
        sel1 = chimera.selection.currentAtoms()
        print(sel1)
        rc("sel #%s:@C6 #%s:@C7 #%s:@C8 #%s:@C9"%(i,i,i,i))
        sel2 = chimera.selection.currentAtoms()
        dihedrals.append((get_dihedral(sel1),get_dihedral(sel2)))
        labels.append(i)

dihedrals = np.array(dihedrals)
x,y = -dihedrals.transpose()
quick_plot(x, y, xlabel=r'$\chi_{1}$', ylabel=r'$\chi_{2}$', labels=labels,
        name="%s/dihedrals.png"%outdir, Type='scatter', fig_size=(12,10),
        xlim=(-180,180), ylim=(-180,180))


#:}}}


# morphing command
#rc("morph movie steps 10")


'''
# If you want to just create images of the structures... This may take longer?

rc('windowsize 600 450')
rc('set depthCue')
rc('set dcColor none')
rc('set dcStart 0.5')
rc('set dcEnd 1.0')
rc('background solid white')
# Apply the publication preset to render a nice image
rc("preset apply publication 1")
#rc('set shadows')
rc('set silhouette')
rc('set silhouetteColor none')
rc('set silhouetteWidth 1.9')



degree = (180./np.pi)*(2./float(turns)*np.pi)
for i in range(0,int(turns)):
    rc('sel #0;focus sel')
    rc('turn y -%s'%degree)
    rc('~sel')
    rc('copy file %s/img_%s.png \
            width %s height %s \
            units %s \
            supersample 4 '%(wd,i,w,h,units))

# FFMPEG
name='img'
num = glob.glob('%s/%s_*.png'%(wd,name))
print num
_list = open('mylist.txt', 'w')
for i in range(0, len(num)):
    run_cmd(" echo 'file %s/%s_%s.png' >> mylist.txt"%(wd,name,i), testing=TESTING)

run_cmd(" echo 'y\n' | ffmpeg -f concat -safe 0 -i %s/mylist.txt -c copy %s/_output.mov "%(wd,wd), testing=TESTING)
run_cmd(" echo 'y\n' | ffmpeg -r %s -i %s/_output.mov -c:v h264 -vf scale=%s:%s %s/output.mov -hide_banner"%(rate,wd,w,h,wd), testing=TESTING)
run_cmd(" echo 'y\n' | ffmpeg -i %s/output.mov %s/output.gif"%(wd,wd), testing=TESTING)

#-vcodec mpeg4


'''







