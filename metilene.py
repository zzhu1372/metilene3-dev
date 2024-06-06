import os
import sys
import time
import argparse
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

metilene_ver = 2.0

###################################################################################################
# Input
###################################################################################################
parser = argparse.ArgumentParser(description='.')
parser.add_argument('-i', "--input",)
parser.add_argument('-o', "--output",)
parser.add_argument('-t', "--threads",)

parser.add_argument('-s', "--skipMetilene",)
parser.add_argument('-om', "--outputImputed",)
parser.add_argument('-u', "--unsupervised",)
parser.add_argument('-gr', "--groupinfo",)
parser.add_argument('-re', "--rerun",)
parser.add_argument('-plt', "--plotTree", default='F')
parser.add_argument('-anno', "--anno",)
parser.add_argument('-gmt', "--gmt",)

parser.add_argument('-m', "--mincpgs", type=int, default=10)
parser.add_argument('-r', "--minDMR", type=int, default=5)
parser.add_argument('-w', "--mindiff", type=float, default=0.1)
parser.add_argument('-wd', "--mindiff_unsup", type=float, default=0.5)
parser.add_argument('-wg', "--mindiff_gsea", type=float, default=0.5)
parser.add_argument('-e', "--mismatch", type=float, default=0.5)

parser.add_argument('-n', "--minN0", type=int, default=1)
parser.add_argument('-g', "--minN", type=int, default=1)
parser.add_argument('-d', "--minNDMR", type=int, default=1)



###################################################################################################
# Install
###################################################################################################
def getMetilene():
    os.system("cd "+os.path.realpath(__file__)[:-len('metilene.py')]+";make")



###################################################################################################
# Run
###################################################################################################
def preprocess(args, headerfile, ifsup, grpinfo=None):
    if ifsup=='unsup':
        cols = pd.read_table(args.input, nrows=0)
        newcols = list(cols.columns)
        # print(newcols)
        for i in range(len(newcols))[2:]:
            newcols[i] = str(i-2)+'_Sample'+str(i-2)
        cols.columns = newcols
        cols.to_csv(headerfile, sep='\t', index=False)
        
    else:
        grp = pd.read_table(grpinfo, index_col='ID')['Group']
        grpid = {}
        j = 0
        for i in sorted(grp.unique()):
            grpid[i] = j
            j += 1
            
        df_grpid = pd.DataFrame(pd.Series(grpid))
        df_grpid.columns = ['Group_ID']
        df_grpid.index.name = 'Group'
        df_grpid.to_csv(args.output+'/'+args.input.split('/')[-1]+'.groupID', sep='\t')
        
        grpdict = grp.map(grpid).to_dict()
        cols = pd.read_table(args.input, nrows=0)
        newcols = list(cols.columns)
        # print(newcols)
        for i in range(len(newcols))[2:]:
            newcols[i] = str(grpdict[newcols[i]])+'_Sample'+str(i-2)#+'_'+newcols[i]
        cols.columns = newcols
        cols.to_csv(headerfile, sep='\t', index=False)


def runMetilene(args, headerfile, ifsup):
    if args.skipMetilene=='T':
        return None
    # print(os.path.realpath(__file__))
    if ifsup=='unsup':
        os.system(os.path.realpath(__file__)[:-3]+" \
                    -t "+str(args.threads)+\
                    " -m "+str(args.mincpgs)+\
                    " -r "+str(args.minDMR)+\
                    " -w "+str(args.mindiff_unsup)+\
                    " -e "+str(args.mismatch)+\
                    " -q "+str(args.mindiff_unsup)+\
                    " -H "+headerfile+\
                    " -d 0 -s 1 -l 1 -p 0 "+\
                    args.input +" > "+\
                    args.output+'/'+args.input.split('/')[-1]+'.unsup.mout' )

    else:
        os.system(os.path.realpath(__file__)[:-3]+" \
                    -t "+str(args.threads)+\
                    " -m "+str(args.mincpgs)+\
                    " -r "+str(args.minDMR)+\
                    " -w "+str(args.mindiff)+\
                    " -e "+str(args.mismatch)+\
                    " -q "+str(args.mindiff)+\
                    " -H "+headerfile+\
                    " -d 0 -s 1 -l 1 -p 0 -O "+str(1*(args.outputImputed=='T'))+' '+\
                    args.input +" > "+\
                    args.output+'/'+args.input.split('/')[-1]+'.aout' )
                    
        if args.outputImputed=='T':
            os.system("grep -v \'//Imputed:\' "+\
            args.output+'/'+args.input.split('/')[-1]+'.aout >' + \
            args.output+'/'+args.input.split('/')[-1]+'.mout')
            
            os.system("head -n1 "+args.input+" > "+\
                        args.output+'/'+args.input.split('/')[-1]+'.imputed')
                        
            os.system("grep \'//Imputed:\' "+\
            args.output+'/'+args.input.split('/')[-1]+".aout|sed \"s/\/\/Imputed://\" >>" + \
            args.output+'/'+args.input.split('/')[-1]+'.imputed')
            
            os.system("rm "+args.output+'/'+args.input.split('/')[-1]+'.aout')
            
        else:
            os.system("mv "+ args.output+'/'+args.input.split('/')[-1]+'.aout ' + \
            args.output+'/'+args.input.split('/')[-1]+'.mout')


def chipseeker(mout, moutPath, anno):
    if anno in ['hg19','HG19']:
        anno = 'TxDb.Hsapiens.UCSC.hg19.knownGene'
    if anno in ['hg38','HG38']:
        anno = 'TxDb.Hsapiens.UCSC.hg38.knownGene'
    cmd = "require("+anno+");require(ChIPseeker);"+\
    "peakfile=\'"+str(os.getcwd())+"/"+moutPath+".bed\';"+\
    "txdb<-"+anno+";"+\
    "peakAnno <- annotatePeak(peakfile, tssRegion=c(-3000, 1000), TxDb=txdb, annoDb=\'org.Hs.eg.db\');"+\
    "write.csv(as.GRanges(peakAnno), \'"+str(os.getcwd())+"/"+moutPath+".bed.csv\')"
    
    mout.sort_values(['chr','start','stop',])[['chr','start','stop']].to_csv(\
    str(os.getcwd())+"/"+moutPath+".bed",sep='\t',index=False,header=None)
    
    os.system('Rscript -e \"'+cmd+'\"')
    
    annoed = pd.read_csv(str(os.getcwd())+"/"+moutPath+".bed.csv", index_col=0)
    
    os.remove(str(os.getcwd())+"/"+moutPath+".bed")
    os.remove(str(os.getcwd())+"/"+moutPath+".bed.csv")
    
    annoed.index = annoed['seqnames']+':'+annoed['start'].astype(str)+'-'+annoed['end'].astype(str)
    annoed['anno'] = annoed['annotation'].apply(lambda x:x.split(' (')[0])
    
    for i in ['distanceToTSS','ENSEMBL','SYMBOL','anno']:
        mout[i] = (mout['chr']+':'+(mout['start']+1).astype(str)+'-'+mout['stop'].astype(str)).map(annoed[i])
        
    return mout


def processOutput(args, ifsup, anno='F'):
    if ifsup=='unsup':
        moutPath = args.output + '/' + args.input.split('/')[-1] + '.unsup.mout'
    else:
        moutPath = args.output + '/' + args.input.split('/')[-1] + '.mout'
    mout = pd.read_table(moutPath)
    mout['meandiffabs'] = mout['meandiff'].apply(abs)

    def rename_cls_pn(x):
        x = x.replace('0','1').replace('4','3')
        if x[0]=='p':
            x = x.replace('1','x').replace('3','1').replace('x','3')
        x = x[1:]
        return x
    mout['sig.comparison'] = ( (1*(mout['meandiff']>0)).map({1:"p", 0:"n"}) \
                                     +mout['sig.comparison']).apply(rename_cls_pn)
    
    mout['#U'] = mout['sig.comparison'].apply(lambda x:(len(x.split('1'))-1))
    mout['#P'] = mout['sig.comparison'].apply(lambda x:(len(x.split('2'))-1))
    mout['#M'] = mout['sig.comparison'].apply(lambda x:(len(x.split('3'))-1))
    
    def calmean(a,b,c):
        a = a.split('|')
        b = b.split('|')
        s = 0
        n = 0
        for i in range(len(a)):
            if b[i]==c:
                s += float(a[i])
                n += 1
        try:
            return s/n
        except:
            return None
            
    mout['meanU'] = mout.apply(lambda x:calmean(x['mean'],x['sig.comparison'],'1'), axis=1)
    mout['meanP'] = mout.apply(lambda x:calmean(x['mean'],x['sig.comparison'],'2'), axis=1)
    mout['meanM'] = mout.apply(lambda x:calmean(x['mean'],x['sig.comparison'],'3'), axis=1)
    # print('# of processed DMRs:',mout.shape[0])
    
    if anno == 'T' and args.anno:
        mout = chipseeker(mout, moutPath, args.anno)
        
    if ifsup=='unsup':
        mout.to_csv(args.output + '/' + args.input.split('/')[-1] + '.unsup.post.mout', \
                    index=False, sep='\t')
    else:
        mout.to_csv(args.output + '/' + args.input.split('/')[-1] + '.post.mout', \
                    index=False, sep='\t')
                    
    return mout



###################################################################################################
# DMR-Freq-based Clustering
###################################################################################################
def recurSplit(arr, ref=0, depth=0, nsep=0, minN=2, minNDMR=100, fulltree=False):
    finalList = []
    depthList = []
    weightList = []
    
    def numVS(a):
        return sorted([a.count('1'),a.count('2'),a.count('3')])[1]
    
    if ref == 0:
        arr = arr.groupby('sig.comparison.bin').sum().sort_values(ascending=False)
        ifSig = 0
        for i in range(len(arr)):
            newref = arr.index[0]
            nsep = arr.iloc[0]
            if nsep > minNDMR and numVS(newref) >= minN:#0.5*nonsep:
                newref = arr.index[i]
                nsep = arr.iloc[i]
                ifSig = 1
                break
        if ifSig:
            finalList.append(newref)
            depthList.append(depth)
            weightList.append(nsep)
        else:
            return None
    else:
        arr = arr.groupby('sig.comparison.bin').sum().sort_values(ascending=False)
        for i in arr.index:
            newref = 0
            if arr[i] < minNDMR:
                return ([],[],[])
            if numVS(i)<minN:
                continue
            newref = i
            finalList.append(newref)
            depthList.append(depth)
            nsep = arr[i]
            weightList.append(nsep)
            break
        # print(newref,fulltree)
        if newref==0 and fulltree:
            # print('FT')
            for i in arr.index:
                newref = 0
                if arr[i] < minNDMR:
                    return ([],[],[])
                if numVS(i)==0:
                    continue
                newref = i
                finalList.append(newref)
                depthList.append(depth)
                nsep = arr[i]
                weightList.append(nsep)
                break
            # print(newref)

    if newref == 0:
        return (finalList, depthList, weightList)
        
    else:
        def mask(x, y, v):
            x = list(x)
            for i in range(len(y)):
                if y[i] != '|':
                    if y[i] != v:
                        x[i] = '0'
            return ''.join(x)
        masked = {}
        for i in ['1','2','3']:
            masked[i] = arr.copy()
            masked[i].index = pd.Series(arr.index).apply(lambda x:mask(x, newref, i))
            resRS = recurSplit(masked[i], mask(newref, newref, i), depth+1, nsep, minN, minNDMR, fulltree)
            finalList += resRS[0]
            depthList += resRS[1]
            weightList+= resRS[2]
            
        return (finalList, depthList, weightList)

    
def plotFTree(cls, reportPath, sids, cmap):
    k = len(cls[0])
    dmrcluster_m = []
    pd.Series(cls[0]).apply(lambda x:dmrcluster_m.append(x.split('|')))
    dmrcluster_m = pd.DataFrame(dmrcluster_m)
    dmrcluster_m = dmrcluster_m.astype(float)
    
    dmrcluster_m.columns = sids
    dmrcluster_m = dmrcluster_m.T
    dmrcluster_m = dmrcluster_m.sort_values(list(range(len(cls[2]))))
    
    ids = list(dmrcluster_m.sort_values(list(range(len(cls[2])))).index+',')+['']
    treestr = ids.copy()
    
    for i in range(len(cls[2])):
        st = {}
        ed = {}
        for j in ['1','2','3']:
            st[j] = dmrcluster_m[i].astype(int).astype(str).sum().find(j)
            ed[j] = dmrcluster_m[i].astype(int).astype(str).sum().rfind(j)
            
            if st[j]!=-1:
                treestr[st[j]] = treestr[st[j]].split(ids[st[j]])[0]+'('+\
                                    ids[st[j]]+treestr[st[j]].split(ids[st[j]])[1]
                treestr[ed[j]] = treestr[ed[j]].split(ids[ed[j]])[0]+ids[ed[j]]+'):'+\
                                    str(cls[2][i])+','+treestr[ed[j]].split(ids[ed[j]])[1]
    
    f = open(reportPath+".nwk", "w")
    f.write(''.join(treestr[:-1]).replace(',)',')')[:-1])
    f.close()
    
    from Bio import Phylo
    import matplotlib.pyplot as plt
    
    tree = Phylo.read(reportPath+".nwk", "newick")

    f,a = plt.subplots(figsize=[20,len(sids)/5])
    Phylo.draw(tree, axes=a, do_show=False, label_colors=cmap)
    plt.savefig(reportPath+'.tree.jpg')
    
    
def plotClustermap(mout, cls, reportPath, sids, finalCls):
    import random
    import seaborn as sns
    import matplotlib.pyplot as plt
    from sklearn.decomposition import PCA
    import numpy as np
    from matplotlib.patches import Patch
    
    k = len(cls[0])
    dmrcluster_m = []
    pd.Series(cls[0]).apply(lambda x:dmrcluster_m.append(x.split('|')))
    dmrcluster_m = pd.DataFrame(dmrcluster_m)
    dmrcluster_m = dmrcluster_m.astype(float)
    
    dmrcluster_m.columns = sids
    dmrcluster_m = dmrcluster_m.T
    dmrcluster_m = dmrcluster_m.sort_values(list(range(len(cls[2]))))
    
    dmrmean_m = [] # for pCA
    mout['mean'].apply(lambda x:dmrmean_m.append(x.split('|')))
    dmrmean_m = pd.DataFrame(dmrmean_m)
    dmrmean_m = dmrmean_m.astype(float).T
    dmrmean_m.index = sids
    pca = PCA(n_components=2)
    X = pd.DataFrame(pca.fit_transform(np.array(dmrmean_m)))
    X.index = dmrmean_m.index

    X['k'] = X[0]-X[0].min()
    X['k'] = X['k']/X['k'].max()
    dmrcluster_m[k] = dmrcluster_m.index.map(X['k'])
    
    def td2(a,b):
        d = 0.01*abs(a[-1]-b[-1])
        for i in range(len(a)-1):
            if a[i]!=b[i]:
                d += max(cls[1])+1 - cls[1][i]
                break
        return d
    
    # print(dmrcluster_m )
    cm = sns.clustermap(dmrcluster_m, metric=td2,\
                        figsize=[3,3],row_cluster=True,col_cluster=False,\
                         dendrogram_ratio=0.2, colors_ratio=0.15, xticklabels=False, yticklabels=False, \
                        method='complete', cmap='Spectral_r')
    lk = cm.dendrogram_row
    # print(lk)
    # print(cls)
    dmrmean_m = []
    
    def editD(a, b):
        a = a.split('|')
        b = b.split('|')
        num = 0
        all = 0
        for i in range(len(a)):
            if a[i] != '0':
                if a[i] != b[i]:
                    num += 1
                all += 1
        if all == 0:
            return 0
        return num/all
        
    denovo_pn2 = mout[['mean','sig.comparison.bin']]
    denovo_pn2['tmp'] = 0
    for i in cls[0]:
        denovo_pn2['tmp'] += 1*(denovo_pn2['sig.comparison.bin'].apply(lambda x:editD(i,x))==0)
    denovo_filtered = denovo_pn2.loc[denovo_pn2['tmp']>0]
    denovo_filtered['mean'].apply(lambda x:dmrmean_m.append(x.split('|')))
    # print(denovo_filtered)
    
    dmrmean_m = pd.DataFrame(dmrmean_m)
    dmrmean_m = dmrmean_m.astype(float)
    dmrmean_m.columns = sids
    dmrmean_m = dmrmean_m[dmrcluster_m.index].T
    
    clsD = finalCls['Group']
    random.seed(0)
    clsCD = {}
    for i in clsD.unique():
        clsCD[i] = (random.randint(1,100)/100,random.randint(1,100)/100,random.randint(1,100)/100)
    cmap = {}
    for i in dmrmean_m.index:
        cmap[i] = clsCD[clsD.to_dict()[i]]

    if dmrmean_m.shape[1]>1:
        cm = sns.clustermap(dmrmean_m,\
            row_colors=[dmrmean_m.index.map(cmap),\
                        (dmrmean_m[0]=='NO').map({False:'white'})],\
            row_linkage=lk.linkage,\
            cmap='Spectral_r', figsize=[0.2*len(sids)+0.1*max([len(i) for i in sids]),0.2*len(sids)], dendrogram_ratio=0.25, xticklabels=False, yticklabels=True, \
            method='ward')
        plt.savefig(reportPath+'.heatmap.jpg')
        plt.savefig(reportPath+'.heatmap.pdf')
    else:
        cm = sns.clustermap(dmrmean_m,\
            row_colors=[dmrmean_m.index.map(cmap),\
                        (dmrmean_m[0]=='NO').map({False:'white'})],\
            row_linkage=lk.linkage,col_cluster=False,\
            cmap='Spectral_r', figsize=[0.2*len(sids)+0.1*max([len(i) for i in sids]),0.2*len(sids)], dendrogram_ratio=0.25, xticklabels=False, yticklabels=True, \
            method='ward')
        plt.savefig(reportPath+'.heatmap.jpg')
        plt.savefig(reportPath+'.heatmap.pdf')

    fig, ax0 = plt.subplots(figsize=(10, 10))
    X['grp'] = X.index.map(cmap)
    # print(X)
    sns.scatterplot(x=X[0],y=X[1],c=X['grp'],ax=ax0,s=30)
    
    label_color_dict = clsCD.copy()
    legend_handles = [Patch(color=color, label=label) for label, color in label_color_dict.items()]
    ax0.legend(handles=legend_handles, ncol=1, )
#    ax0.axis('off')
    plt.savefig(reportPath+'.dmrpca.jpg')

    return cmap


def clustering(mout, args):
    minN0 = args.minN0
    minN = args.minN
    minNDMR = args.minNDMR
    mindiff_unsup = args.mindiff_unsup

    def rename_cls_pn2(x):
        if x.count('3')>x.count('1'):
            x = x.replace('1','2')
        elif x.count('3')<x.count('1'):
            x = x.replace('3','2')
        return x
    mout = mout.loc[(mout['#U']>=minN0)\
                    &(mout['#M']>=minN0)\
                    &(mout['meandiffabs']>mindiff_unsup)]
    mout['sig.comparison.bin'] = mout['sig.comparison'].apply(rename_cls_pn2)
    ranked = mout[['sig.comparison.bin','meandiffabs']].groupby('sig.comparison.bin').\
        sum()['meandiffabs'].sort_values(ascending=False)
    cls = recurSplit(ranked.sort_values(ascending = False), minN=minN, minNDMR=minNDMR)
    # print(cls)
    if cls is None:
        return None

    reportPath = args.output+'/'+args.input.split('/')[-1]
    sids = list(pd.read_table(args.input, nrows=0).columns[2:])
    
    finalCls = pd.DataFrame([i.split('|') for i in cls[0]]).sum()
    finalCls.index = sids
    rename_cls_id = {}
    j=0
    for i in sorted(finalCls.unique()):
        rename_cls_id[i] = 'G'+str(j)
        j+=1
    finalCls = finalCls.map(rename_cls_id)
    finalCls = pd.DataFrame(finalCls)
    finalCls.columns = ['Group']
    finalCls.index.name = 'ID'
    finalCls.to_csv(reportPath+'.clusters', sep='\t')
    
    cmap = plotClustermap(mout, cls, reportPath, sids, finalCls)
    if args.plotTree=='T':
        cls_full = recurSplit(ranked.sort_values(ascending = False), \
                              minN=minN, minNDMR=0, fulltree=True)
        plotFTree(cls_full, reportPath, sids, cmap)
    
    return finalCls



###################################################################################################
# Run
###################################################################################################
def report(args, start_time, end_time, unmout, finalCls, mout):
    with open(os.path.realpath(__file__)[:-len('metilene.py')]+'template.html', 'r') as template_file:
        template_content = template_file.read()

    final_html = template_content.replace('<h2>Metilene Report for XXX</h2>', '<h2>Metilene Report for '+args.input.split('/')[-1]+'</h2>')
    final_html = final_html.replace('<div>Version: XXX</div><br>', '<div>Version: '+str(metilene_ver)+'</div><br>')
    final_html = final_html.replace('<div>Command: XXX</div><br>', '<div>Command: '+''.join([i+' ' for i in sys.argv])+'</div><br>')
    final_html = final_html.replace('<div>Parameters: XXX</div><br>', '<div>Parameters: <br>'+str(args).split('Namespace')[-1][1:-1]+'</div><br>')
    final_html = final_html.replace('<div>Start time: XXX</div><br>', '<div>Start time: '+str(start_time)+'</div>')
    final_html = final_html.replace('<div>End time: XXX</div><br>', 'End time: '+str(end_time)+'</div><br>')

    final_html = final_html.replace('<div>Number of unsupervised DMRs: XXX</div><br>', 'Number of unsupervised DMRs: '+str(unmout.shape[0])+'</div><br>')
    
    final_html = final_html.replace('<div>Number of clusters: XXX</div><br>', 'Number of clusters: '+str(len(finalCls['Group'].unique()))+'</div><br>')
    cls_table_html = finalCls.to_html(escape=False)
    final_html = final_html.replace('<div id="pandas_table_placeholder_cluster"></div>', cls_table_html)

    if args.plotTree=='T':
        final_html = final_html.replace('./clstree.png', './'+args.input.split('/')[-1]+'.tree.jpg')
    else:
        final_html = final_html.replace('./clstree.png', '')
        
    final_html = final_html.replace('./heatmap.png', './'+args.input.split('/')[-1]+'.heatmap.jpg')
    final_html = final_html.replace('./dmrpca.png', './'+args.input.split('/')[-1]+'.dmrpca.jpg')

    final_html = final_html.replace('<div>Number of supervised DMRs: XXX</div><br>', 'Number of supervised DMRs: '+str(mout.shape[0])+'</div><br>')

    table = mout.groupby('sig.comparison').count()[['p']].sort_values('p', ascending=False)[:10]
    table.columns = ['#DMRs']
    
    def decodeSigCmp(x):
        upmlist = {'1':[], '2':[], '3':[]}
        x = x.split('|')
        grp_dict = pd.read_table(args.output+'/'+args.input.split('/')[-1]+'.groupID',\
                                  index_col='Group_ID')
        grp_dict = grp_dict['Group'].to_dict()
        for i in range(len(x)):
            upmlist[x[i]].append(grp_dict[i])
        return 'Unmet:'+','.join(upmlist['1'])+' - vs - '+'Met:'+','.join(upmlist['3'])
    
    gseapopup = ''
    if args.gmt and args.anno:
        import gseapy as gp
        j = 0
        for i in table.index:
            gene_sets = args.gmt
            gene_list = list(set(mout.loc[(mout['sig.comparison']==i)&(mout['meandiffabs']>args.mindiff_gsea)]['SYMBOL'].dropna()))
            
            try:
                enr = gp.enrichr(gene_list=gene_list,
                            gene_sets=gene_sets,
                            organism='human',
                            outdir=args.output+'/'+args.input.split('/')[-1]+'.gsea/'+i.replace('|','_'),
                            cutoff = 1,
                            format = 'jpg',
                            )
            except:
                print("GSEA error:",gene_list)
                        
            fig_path = './'+args.input.split('/')[-1]+'.gsea/'+i.replace('|','_')+\
                        "/"+args.gmt.split('/')[-1]+".human.enrichr.reports.jpg"
                        
            gseapopup += "<div id=\"popupgsea"+str(j)+"\" class=\"popup\"><br>\
                <button onclick=\"hidePopup('popupgsea"+str(j)+"')\">Close</button>\
                    <p>GSEA for "+decodeSigCmp(i)+"</p><img src="+fig_path+" height=\"200\"><br></div>\n"
            j+=1
            
        table['GSEA'] = ["<button onclick=\"showPopup('popupgsea"+str(i)+\
                        "')\">Click to show GSEA results</button>" \
                        for i in range(table.shape[0])]
    
    table.index = [decodeSigCmp(i) for i in table.index]

    table_html = table.to_html(escape=False)
    final_html = final_html.replace('<div id="pandas_table_placeholder_dmr"></div>', table_html)
    final_html = final_html.replace('<div id="gsea_placeholder"></div>', gseapopup)
    
    with open(args.output+'/'+args.input.split('/')[-1]+'.report.html', 'w') as final_file:
        final_file.write(final_html)



###################################################################################################
# main
###################################################################################################
def main():
    start_time = time.ctime()
    
    args = parser.parse_args()
    print(args)
    
    getMetilene()
    try:
        os.mkdir(args.output)
    except:
        pass
        
    if args.unsupervised=='T':
        headerfile = args.output+'/'+args.input.split('/')[-1]+'.unsup.header'
        preprocess(args, headerfile, 'unsup')
        runMetilene(args, headerfile, 'unsup')
        unmout = processOutput(args, 'unsup')
        finalCls = clustering(unmout, args)
        if finalCls is None:
            print('ERROR: No cluster found. Please check the data or use smaller meandiff for clustering.')
            return None

        if args.rerun=='F':
            return None
        
        headerfile = args.output+'/'+args.input.split('/')[-1]+'.header'
        preprocess(args, headerfile, 'sup', \
                   args.output+'/'+args.input.split('/')[-1]+'.clusters')
        runMetilene(args, headerfile, 'sup')
        mout = processOutput(args, 'sup', anno='T')

        end_time = time.ctime()
        report(args, start_time, end_time, unmout, finalCls, mout)
    
    else:
        headerfile = args.output+'/'+args.input.split('/')[-1]+'.header'
        preprocess(args, headerfile, 'sup', \
                   args.groupinfo)
        runMetilene(args, headerfile, 'sup')
        mout = processOutput(args, 'sup', anno='T')
        
        
main()
