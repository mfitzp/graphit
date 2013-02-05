from numpy.random import uniform, seed
from collections import OrderedDict
from matplotlib.mlab import griddata
from matplotlib import ticker
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import re,os
import random
import csv

from optparse import OptionParser
from collections import defaultdict

import utils

parser = OptionParser()

parser.add_option("-m", "--metabolite", dest="metabolite", default='',
                  help="name of metabolite to plot graph for")

parser.add_option("-f", "--file", dest="file", default='plsdump.csv',
                  help="csv file containing dumped data from plsmulti", metavar="FILE")

parser.add_option("-b", "--batch", action="store_true", dest="batch_mode", default=False,
                  help="batch mode, process all metabolite matches with same parameters")

parser.add_option("-l", "--longitude", dest="longitudinal", default='(\d+)',
                  help="regex pattern for longitudinal data: grouped area is x axis")

parser.add_option("-c", "--control", dest="control", default=False,
                  help="regex pattern for control data")

parser.add_option("--multiplot", action="store_true", dest="multiplot", default=False,
                  help="plot more than one graph on a figure")

parser.add_option("--shareyaxis", action="store_true", dest="shareyaxis", default=False,
                  help="whether to share y axis scale (and show on multiplots)")

parser.add_option("-s", "--search", dest="search", default=None,
                  help="only show classes matching this regex")

parser.add_option("--styles", dest="styles", default=None,
                  help="style pattern to use to assign colors/markers/lines to lines; (3),(regex),(groups)")

parser.add_option("-d", "--display", action="store_true", dest="display", default=False,
                  help="display resulting graphs to screen")

parser.add_option("--xlabel", dest="xlabel", default='',
                  help="x axis label for graph")

parser.add_option("--ylabel", dest="ylabel", default='',
                  help="y axis label for graph")

parser.add_option("--logx", dest="logx", action="store_true", default=False,
                  help="log x axis")

parser.add_option("--logy", dest="logy", action="store_true", default=False,
                  help="log y axis")

parser.add_option("--smooth", dest="smooth", action="store_true", default=False,
                  help="smooth graph through cubic interpolation")
                  
parser.add_option("--ylim", dest="ylim", default=None,
                  help="min,max for y axis")
                  
parser.add_option("-t", "--title", dest="title", default=None,
                  help="title for graph")

parser.add_option("--dpi", dest="dpi", default=72,
                  help="dpi for output")

parser.add_option("--format", dest="format", default='png',
                  help="fileformat for output")

parser.add_option("--annotate", action="store_true", dest="annotate", default=False,
                  help="show command annotation for generation")


(options, args) = parser.parse_args()


colors = ['r','b','g','c','m','y','k']
markers = ['o', 's','v','^','D','+','x']
linestyles = ['solid', 'dashed', 'dashdot'] # Do not include dotted, used for controls

if options.styles is None:
    # Build a full table we'll just iterate over
    styles = [(z,y,x) for x in linestyles for y in markers for z in colors]

# Extract file root from the edge file name
filebase = os.path.splitext(options.file)[0]
[null, sep, string ] = filebase.partition('-')
filesuff = sep + string
nodes = OrderedDict()

(metabolites, allquants, globalylim) = utils.read_metabolite_datafile(options.file, options)

# Turn off interactive plotting, speed up
plt.ioff()
figures = list()
multiax = None
ymaxf = 0

for metabolite in metabolites[:]:
    print "Processing %s" % metabolite
    quants = allquants[metabolite]

    if options.search:
        okeys = quants.keys()
        for label in quants.keys():
            match = re.search(options.search, label) 
            if not match:
                del quants[label]
        print "Filter matching classes '%s' with '%s' gives '%s'" % (', '.join(okeys), options.search, ', '.join( quants.keys() ) )
        if len( quants.keys() ) == 0:
            print "Nothing left!; deleting metabolite"
            metabolites.remove( metabolite )
            continue
    
    
    # Apply regex to class string for each variable;
    # the *.(\d+)*. longitudinal bit extracts that
    # use remainder for class assignment *at that timepoint*
    # Possible to use the same structure?
    
    timepoints = set()
    classes = set()
    classtrans = defaultdict(list)
    non_decimal = re.compile(r'[^\d.]+')
    
    for label in quants.keys():
        match = re.search(options.longitudinal, label)
        if match:
            # Rewrite the class label to not include the match content  
            classlabel = label.replace(match.group(1),'')
            timepoint = float( non_decimal.sub('', match.group(1)) )
            timepoints.add( timepoint ) # Build timepoint list
            classes.add( classlabel ) # Build class list
            
            classtrans[ (classlabel, timepoint) ] = label # Store translation for lookup
            
    if len(timepoints) ==0:
        print "No matching classes found for longitude regex, try again"
        exit()

    ind = list()
    graphs = dict()
    classes = list(classes)

    # Remove axis timepoint duplicates and sort
    timepoints = list( timepoints )
    timepoints.sort(key=float)
    print "Axis longitude completed: " + ", ".join(str(x) for x in timepoints)
    print "Re-classification: " + ", ".join(classes)

    
    # If set run styles setting against the class list and build a lookup table of variants
    # with assigned colours and styles
    if options.styles:
        
        classstyles = dict()
        classmatch = defaultdict(list)
        stylesets = list()
        for n, styleret in enumerate( options.styles.split(',') ):
            stylere = re.compile('(%s)' % styleret)
            stylesets.append( set() )
            # Build table, then assign
            for classl in classes:
                match = stylere.search(classl)
                if match:
                    stylesets[n].add(match.group(1))
                    classmatch[classl].append( match.group(1) )
                else:
                    classmatch[classl].append( None )
    
        # Now have 3 sets of styles, assign
        for (classl,classm) in classmatch.items():
            classstyles[classl] = (
                colors[ list( stylesets[0] ).index(classm[0]) ], 
                markers[ list( stylesets[1] ).index(classm[1]) ], 
                linestyles[ list( stylesets[2] ).index(classm[2]) ],
                )
                
    # Get common substring from classes, for improved naming of saved graph file
    common_classname = ''
    # Exclude classes matching control regexp
    # ...
    classnc = [classl for classl in classes if not re.search(".*%s.*" % options.control, classl)]
    
    if len(classnc) > 1 and len(classnc[0]) > 0:
        for i in range(len(classnc[0])):
            for j in range(len(classnc[0])-i+1):
                if j > len(common_classname) and all(classnc[0][i:i+j] in x for x in classnc):
                    common_classname = classnc[0][i:i+j]
    common_classname = common_classname.strip("-")


    # Output is a dict of lists, containing values matching each point on the X axis (or None if not existing)
    ymin = 0
    ymax = 0
    controls = 0
    
    datasets = defaultdict(list)

    for n, classl in enumerate(classes):
        graph = {
            'timepoints':list(),
            'means':list(),
            'stddev':list(),
            'samples':list(),
            'color':'k',
            'control':False,
        }
        
        if options.control:
               match = re.search(".*%s.*" % options.control, classl)
               if match: # This is a control sample, mark as such
                    graph['control'] = True
                    controls += 1
        
        for t in timepoints:
            
            if (classl,t) in classtrans:
                key = classtrans[ (classl, t) ] # Get original classname for this timepoint/class combo
                graph['timepoints'].append( float(t) )
                graph['means'].append( np.mean( quants[key]) )
                graph['stddev'].append( np.std( quants[key]) )
                graph['samples'].append( len(quants[key] ) )
                ymax = max( ymax, max( quants[key] ) )
                ymin = min( ymin, min( quants[key] ) )
        
        if options.styles:
            graph['style'] = classstyles[classl]
        else:
            graph['style'] = styles[n]
        
        # Store completed graph
        graphs[classl] = graph


    if options.shareyaxis:
        ylim = globalylim
    else:
        ylim = (ymin, ymax)
    
    
    if options.multiplot:
        # Keep using same figure, append subplots
        if not multiax:
            #adjustprops = dict(left=0.1, bottom=0.1, right=0.97, top=0.93, wspace=0.2, hspace=0.2) 
            fig = plt.figure()
            multiax = fig.add_subplot(1, len(metabolites), 1 )
            figures.append( multiax ) 
        else:
            if not options.multiaxes:
                fp = fig.add_subplot(1, len(metabolites), len(figures)+1, sharey=multiax )
                plt.setp(fp.get_yticklabels(), visible=False)
                plt.setp(fp.get_yaxis(), visible=False)
            else:
                fp = fig.add_subplot(1, len(metabolites), len(figures)+1 )

            figures.append( fp ) 
    else:
        # New figure
        figures.append( plt.figure() )

    for classl in classes:
        graph = graphs[classl]
        
        if options.smooth and len(graph['timepoints']) >= 3:    

            l = plt.errorbar( graph['timepoints'], graph['means'], yerr=graph['stddev'], label=classl, fmt=graph['style'][0], marker=graph['style'][1], linestyle=graph['style'][2] )
            # Calculate spline and plot that over
            # Duplicate data to additional timepoints
            # Convert data to linear timecourse first (doh!)
            ynew = list()
            for n in range(len(graph['means'][:-1])):
                t1 = graph['timepoints'][n]
                t2 = graph['timepoints'][n+1]
                m1 = graph['means'][n]
                m2 = graph['means'][n+1]
                dm = int(t2-t1)
                for tp in range(dm):
                    av = ( (m1*(dm-tp)/dm) + (m2*tp)/dm )
                    ynew.append( av )
                    
            window_len = 11
            
            xnew = np.linspace( min(graph['timepoints']), max(graph['timepoints']), len(ynew) )
            ynew = np.array( ynew )
            
            # Iteratively smooth with a small window to minimise distortion
            # smooth the xaxis alongside to prevent x-distortion
            for x in range( 20 ):
                ynew = utils.smooth(ynew, window_len=window_len, window='blackman')
                xnew = utils.smooth(xnew, window_len=window_len, window='blackman')
            
            plt.plot( xnew, ynew, color=graph['color'])
        else:
            if graph['control']:
                plt.fill_between( graph['timepoints'], [a+b for a, b in zip( graph['means'], graph['stddev'] )], y2=[a-b for a, b in zip( graph['means'], graph['stddev'] )], alpha=0.05, color=graph['style'][0])
                plt.plot( graph['timepoints'], graph['means'], color=graph['style'][0], linestyle=graph['style'][2], label=classl, alpha=0.2)
            else:       
                plt.errorbar( graph['timepoints'], graph['means'], yerr=graph['stddev'], label=classl, fmt=graph['style'][0], marker=graph['style'][1], linestyle=graph['style'][2])


    if options.control:
        pass                
        #for classl in baselines:
        #    p = mpatches.Rectangle(xy, width, height, facecolor="orange", edgecolor="red")
        #    plt.gca().add_patch(p)
        #    plt.draw()       
    
    ax = plt.gca()
    handles, labels = ax.get_legend_handles_labels()
    hl = sorted( zip(labels,handles) )
    labels =  [hi[0] for hi in hl]
    handles = [hi[1] for hi in hl]
    ax.legend(handles, labels)   
    
    # Optional log x and y axes
    if options.logx:
        ax.set_xscale('log')

    if options.logy:
        ax.set_yscale('log')


    #plt.xticks(timepoints, timepoints, rotation=45)

    if options.title:
        plt.title("%s (%s)" % (options.title,metabolite))
    else:
        plt.title(metabolite)
    
    plt.gca().xaxis.set_label_text(options.xlabel)
    plt.gca().yaxis.set_label_text(options.ylabel)
    
    # Add some padding either side of graphs
    #plt.xlim( ind[0]-1, ind[-1]+1)
    
    if options.annotate:
        utils.annotate_plot(plt,options)
    
    fig = plt.gcf()
    # Horizontal axis through zero
    #plt.axhline(0, color='k')

    if options.multiplot: #Scale the multiplots a bit more reasonably
        fig.set_size_inches(5 + len(metabolites)*3,6)
    else:
        fig.set_size_inches(8,6)
    
    # Get ylimits for significance bars

    if options.ylim:
        ylim=options.ylim.split(',')
        plt.ylim( int(ylim[0]), int(ylim[1]) )
        
    if not options.multiplot:   
        # Adjust plot on multiplot
        #plt.ylim( ymin, ymaxf)  
        ymaxf = 0
        print "Save as 'long%s-%s-%s.%s'" % (filesuff, metabolite, common_classname, options.format)
        plt.savefig('long%s-%s-%s.%s' % (filesuff, metabolite, common_classname, options.format), dpi=options.dpi, transparent=False )

if options.multiplot:   
    # Adjust plot on multiplot
    plt.ylim( ymin, ymaxf )  
    print "Save as 'long%s-%s-%s.%s'" % (filesuff, '-'.join(metabolites), common_classname, options.format)
    plt.savefig('long%s-%s-%s.%s' % (filesuff, '-'.join(metabolites), common_classname, options.format), dpi=options.dpi, transparent=False )

if options.display:
    plt.show()   
plt.close()