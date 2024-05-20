import os
import argparse
import pandas as pd


###################################################################################################
# Input
###################################################################################################
parser = argparse.ArgumentParser(description='.')
parser.add_argument('-i', "--input",)
parser.add_argument('-o', "--output",)
parser.add_argument('-t', "--threads",)

parser.add_argument('-s', "--skipMetilene",)

parser.add_argument('-m', "--mincpgs",)
parser.add_argument('-r', "--minDMR",)
parser.add_argument('-w', "--mindiff",)
parser.add_argument('-e', "--minDMR2",)
parser.add_argument('-q', "--mindiff2",)

parser.add_argument('-n', "--minN0", type=int)
parser.add_argument('-g', "--minN", type=int)
parser.add_argument('-d', "--minNDMR", type=int)


###################################################################################################
# Install
###################################################################################################
def getMetilene():
    pass

###################################################################################################
# Run
###################################################################################################
def runMetilene(args):
    print(args.skipMetilene)
    if (args.skipMetilene=='F'):
        os.system("cd ./; \
                    ./metilene \
                    -t "+str(args.threads)+\
                    " -m "+str(args.mincpgs)+\
                    " -r "+str(args.minDMR)+\
                    " -w "+str(args.mindiff)+\
                    " -e "+str(args.minDMR2)+\
                    " -q "+str(args.mindiff2)+\
                    " -d 0 -s 1 -l 1 -p 1 "+\
                    args.input +" > "+\
                    args.output+'/'+args.input.split('/')[-1]+'.mout' )

def processOutput(args):
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
    print(mout.shape)
    mout.to_csv(args.output + '/' + args.input.split('/')[-1] + '.post.mout')
    return mout


###################################################################################################
# DMR-Freq-based Clustering
###################################################################################################
def recurSplit(arr, ref=0, depth=0, nsep=0, minN=2, minNDMR=100):
    finalList = []
    depthList = []
    weightList = []
    
    def numVS(a):
        return sorted([a.count('1'),a.count('2'),a.count('3')])[1]
    
    if ref == 0:
        arr = arr.groupby('sig.comparison.bin').sum().sort_values(ascending=False)
        for i in range(len(arr)):
            newref = arr.index[i]
            nsep = arr.iloc[i]
            if nsep > minNDMR and numVS(newref) >= minN:#0.5*nonsep:
                break
        finalList.append(newref)
        depthList.append(depth)
        weightList.append(nsep)
        
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
            resRS = recurSplit(masked[i], mask(newref, newref, i), depth+1, nsep, minN, minNDMR)
            finalList += resRS[0]
            depthList += resRS[1]
            weightList+= resRS[2]
            
        return (finalList, depthList, weightList)
    
def plotFTree(cls, reportPath, sids):
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
                treestr[st[j]] = treestr[st[j]].split(ids[st[j]])[0]+'('+ids[st[j]]+treestr[st[j]].split(ids[st[j]])[1]
                treestr[ed[j]] = treestr[ed[j]].split(ids[ed[j]])[0]+ids[ed[j]]+'):'+str(cls[2][i])+','+treestr[ed[j]].split(ids[ed[j]])[1]
    
    f = open("test.nwk", "w")
    f.write(''.join(treestr[:-1]).replace(',)',')')[:-1])
    f.close()
    
    from Bio import Phylo
    import matplotlib.pyplot as plt
    
    tree = Phylo.read("test.nwk", "newick")

    f,a = plt.subplots(figsize=[30,15])
    Phylo.draw(tree, axes=a)
    plt.savefig(reportPath+'/tree.jpg')

def clustering(mout, args):
    minN0 = args.minN0
    minN = args.minN
    minNDMR = args.minNDMR

    def rename_cls_pn2(x):
        if x.count('3')>x.count('1'):
            x = x.replace('1','2')
        elif x.count('3')<x.count('1'):
            x = x.replace('3','2')
        return x
    mout = mout.loc[(mout['#U']>=minN0)\
                    &(mout['#M']>=minN0)\
                    &(mout['meandiffabs']>0.5)]
    mout['sig.comparison.bin'] = mout['sig.comparison'].apply(rename_cls_pn2)
    ranked = mout[['sig.comparison.bin','meandiffabs']].groupby('sig.comparison.bin').\
        sum()['meandiffabs'].sort_values(ascending=False)
    cls = recurSplit(ranked.sort_values(ascending = False), minN=minN, minNDMR=minNDMR)
    print(cls)

    reportPath = args.output
    sids = list(pd.read_table(args.input, nrows=0).columns[2:])
    plotFTree(cls, reportPath, sids)

    return cls


        
###################################################################################################
# main
###################################################################################################
def main():
    args = parser.parse_args()

    getMetilene()

    runMetilene(args)
    
    mout = processOutput(args)

    clustering(mout, args)
    
main()